

function [RealTimeTraversing, flag] = identify_single_v2_6(MarkerSetsNK, mars, dism1m2)

% 识别
% 输入单模板的边信息、点结构体以及捕捉（采集）点的边信息，返回识别成功与否(flag)和模板点与实时点之间的识别匹配关系(RealTimeTraversing)

% identify_v1_3         使用.mars文件中的排序方式，按模板边的顺序遍历识别，识别对象仍为模板点
% identify_v1_4         添加模板中大于1个连通图的情况的处理
% identify_v1_5         将识别过程中的栈 TraversalStack 设为全局变量
% identify_v1_6         整理整个识别过程，明确变量名，去掉不需要的结构体
% identify_v2_0         非刚体、多模板、单目标识别，模板及采集坐标信息从外部输入
% identify_single_v2_1  改名为identify_single，单模板单目标识别。多模板多目标由identify_whole实现。
% identify_single_v2_2  修改cpt_nextRTT中寻找当前实时点连接的其他实时点的方法，优化出栈入栈的过程。
% identify_single_v2_3  在步骤“确定公共点为当前模板点的备选匹配点后，为另一个端点的模板点提供备选点”的同时更新模板边的备选边
% identify_single_v2_4  调整“确定公共点为当前模板点的备选匹配点后，为另一个端点的模板点提供备选点”的同时更新模板边的备选边的方法
%                       更新模板边放在更新端点备选点一起，删除变量 RealtimeSets.
% identify_single_v2_5  将按点归类步骤中，更新模板点和模板边的内容分开，更新模板边的步骤放到全部计算完模板点的备选点之后。
% identify_single_v2_6  调整按模板边归类时的跳出条件，修改注释与识别文档中的名称对应。

flag = 1;  % 1,识别成功；0,识别失败
RealTimeTraversing = [] ;       % 最后一个元素表示当前正在识别的模板边的编号，
                                % 倒数第二个元素表示当前正在识别的模板点编号，
                                % 从第一个到倒数第二个元素表示模板点对应的实时点号，
                                % 如 RealTimeTraversing(3)=48 表示3号模板点对应的实时点号为48号

% tic

NumberOfSide = mars.NumberOfSide ;
NumberOfMarker = mars.NumberOfMarker ;
ExtraStretch = mars.ExtraStretch ; %误差浮动范围
TemplatePoint = mars.TemplatePoint ; %模板点的结构体

% -----------------------  2.2 按边长归类  --------------------------------
% 将采集边按模板边的长度归类，筛选模板边的备选匹配边
Match_Side = cell(NumberOfSide,1) ; % Match_TempLateSide_RealTimeSide 第i个元素存储与第i号模板边匹配的采集边

i=1; % index_MarkerSetsNK
j=1; % index_dism1m2
Nj = size(dism1m2,1) ;
MarkerSetsNK = sortrows(MarkerSetsNK) ;
dism1m2 = sortrows(dism1m2) ;
while( i<=NumberOfSide  ) 
    SideNum = MarkerSetsNK(i,5) ; %模板边的编号
    
    % 寻找与当前模板边长度匹配的采集边
    if dism1m2(j,1) < MarkerSetsNK(i,1) - ExtraStretch 
        j = j+1 ;
    else
        tempMatch = [] ;
        tj = j ;
        while (dism1m2(tj,1) < MarkerSetsNK(i,2) + ExtraStretch) && tj<Nj
            tempMatch = [tempMatch; dism1m2(tj,:)] ;
            tj = tj+1 ;
        end
        
        if isempty(tempMatch) || j>Nj
%             fprintf('%s\n',['构造边的匹配关系的结构体，模板template中第',int2str(SideNum),'号边没有匹配到~']) ;
            flag = 0 ;
            break ;
        end
        
        Match_Side{SideNum} = tempMatch ;
        i = i+1 ;
    end
end
if flag<1, return; end

% t1 = toc ;
% fprintf('按长度筛选备选边：%f \n', t1) ;

% ---------------------------  2.2 end  ----------------------------------

% ------------------------  2.3 按点归类  ---------------------------------
% 匹配一下试试看呢…… = =
% 一个一个地计算模板点的备选点~
Match_Point = cell(NumberOfMarker,1) ; % Match_TempLatePoint_RealTimePoint 模板点与实时点的匹配关系. 第i行存放与第i号模板点匹配的实时点编号

for i=1:NumberOfMarker
    MarkerNum = TemplatePoint(i).Num ;
    LinkageNum = TemplatePoint(i).LinkageNum ;
    if isempty(LinkageNum), continue; end
    
    % Step2.3.1 同一个模板点连接的各模板边的备选边必定含有一个公共端点
    if isempty(Match_Point{i})
        p_re = Match_Side{LinkageNum(1)}(:,2:3) ; %p_re, 与模板点连接的备选边上存在的公共点
        p_re = unique(p_re(:)) ;
        start = 2 ;
    else
        p_re = Match_Point{i} ;
        start = 1 ;
    end
    
    for j = start:length(LinkageNum)
        temp = Match_Side{LinkageNum(j)}(:,2:3) ;
        p_re = intersect(p_re(:),temp(:)) ;
    end
    
    if isempty(p_re) %若不存在公共点，则说明此模板点没有匹配的实时点，此帧不能整体识别
%         fprintf('%s\n',['构造模板点与实时点的匹配关系，模板template中第',int2str(MarkerNum),'号点没有匹配到~']) ;
        flag = 0 ;
        break ;
    end
    
    Match_Point{i} = p_re ;

    % Step2.3.2 确定公共点为当前模板点的备选匹配点后，为另一个端点的模板点提供备选点
    LinkagePoint = TemplatePoint(i).LinkagePoint ;
    for ilp = 1:length(LinkagePoint) % index_LinkagePoint 按连接边的顺序依次比较
        
        % 寻找与该连接边另一端点匹配的实时点和对应的采集边
        tms = Match_Side{LinkageNum(ilp)} ; %temp_Match_Side
        nowMP = [] ; %new_Match_Point
        for ipre = 1:length(p_re)
            tpre = p_re(ipre) ;
            temp = tms(tms(:,2)==tpre | tms(:,3)==tpre, :) ;
            
            temp = temp(:,2:3) ;
            temp(temp==tpre) = [] ;
            nowMP = [nowMP; temp(:)] ; 
        end
        nowMP = unique(nowMP) ;
        
        % 更新点和边的匹配信息
        % 将新算出来的点、边与匹配表中的点、边匹配信息比较，留下两边都存在的项目
        if isempty(Match_Point{LinkagePoint(ilp)})
            Match_Point{LinkagePoint(ilp)} = nowMP ;
        else
            newMP = intersect(Match_Point{LinkagePoint(ilp)},nowMP) ;   %new_Match_Point
            if isempty(newMP)
%                 fprintf('%s\n',['i=',int2str(MarkerNum),',模板template中第',int2str(LinkagePoint(ilp)),'号点没有匹配到~']) ;
                flag = 0 ;
                break ;
            else
                Match_Point{LinkagePoint(ilp)} = newMP ;
            end
        end
    end
    
    if flag<1, break; end
end %for i=1:NumberOfMarker
if flag<1, return; end

% ---------------------------  2.3 end  ----------------------------------

% --------------------  2.4 更新模板边对应的备选边  ----------------------------

MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
for i = 1:NumberOfSide
    TemplateMarker1 = MarkerSetsNK(i,3) ;
    TemplateMarker2 = MarkerSetsNK(i,4) ;
    RealtimeMarker1 = Match_Point{TemplateMarker1} ;
    RealtimeMarker2 = Match_Point{TemplateMarker2} ;
    
    tempMatch = Match_Side{i} ;
    newMatch = [] ;
    for j = 1:length(RealtimeMarker1)
        for k = 1:length(RealtimeMarker2)
            idx = tempMatch(:,3)==RealtimeMarker1(j) & tempMatch(:,2)==RealtimeMarker2(k) ;
            newMatch = [newMatch; tempMatch(idx, :) ];
            idx = tempMatch(:,2)==RealtimeMarker1(j) & tempMatch(:,3)==RealtimeMarker2(k) ;
            newMatch = [newMatch; tempMatch(idx, :) ];
        end
    end
    
    if isempty(newMatch)
%         fprintf('%s\n',['为模板边更新备选边，模板template中第',int2str(i),'号边没有匹配到~']) ;
        flag = 0 ;
        break ;
    else
        Match_Side{i} = newMatch ;
    end
end
if flag<1, return; end

% ---------------------------  2.4 end  ----------------------------------

% --------------------  2.5 按模板边顺序遍历识别  --------------------------

% 设个存档点试试看呢……
% 用模板边来匹配~？

MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
idxTLpoint = NumberOfMarker+1 ; % RealTimeTraversing 中存放模板点编号的位置
idxTLside  = NumberOfMarker+2 ; % RealTimeTraversing 中存放模板边编号的位置
LengthRTT  = idxTLside ;      % RealTimeTraversing 的长度
RealTimeTraversing = zeros(1,LengthRTT) ;       % 最后一个元素表示当前正在识别的模板边的编号，
                                                % 倒数第二个元素表示当前正在识别的模板点编号，
                                                % 从第一个到倒数第二个元素表示模板点对应的实时点号，
                                                % 如 RealTimeTraversing(3)=48 表示3号模板点对应的实时点号为48号
TraversalStack = [] ;

% 2.5.1 初始化遍历向量
% 模板边连接的编号最小的模板点作为第一个模板点，即识别的起点，计算第一个模板点对应的实时点编号
[RealTimeTraversing, TraversalStack] = InitRTpoint( 0, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, Match_Point) ;

% 2.5.2 计算下一个模板点的识别点
% 按模板边的排列顺序，依次计算后续模板点对应的实时点编号
idxMarkerSets = 1 ;
while idxMarkerSets <= NumberOfSide
    
    % 计算下一个可能出现的RealTimeTraversing
    [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, MarkerSetsNK, idxTLpoint, idxTLside, ExtraStretch, Match_Point, Match_Side);
    
    nRTT = size(nextRTT,1) ; %nextRTT中含有的过程数
    if nRTT == 1
        RealTimeTraversing = nextRTT ;
    elseif nRTT < 1
        n = size(TraversalStack,1) ;
        if n < 1
%             disp('整体识别失败~') ;
            flag = 0 ;
            break ;
        end
        % 出栈，取栈顶的过程作为当前过程
        RealTimeTraversing = TraversalStack(n ,:) ;
        TraversalStack(n ,:) = [] ;     
    else
        RealTimeTraversing = nextRTT(nRTT,:) ;
        % 入栈
        TraversalStack = [TraversalStack; nextRTT(1:nRTT-1,:)] ;
    end
    idxMarkerSets = RealTimeTraversing(idxTLside) ;
end
if flag<1, return; end

RealTimeTraversing = RealTimeTraversing' ;

% t2 = toc ;
% fprintf('识别：%f \n', t2-t1) ;
% t1 = t2 ;

% ---------------------------  2.5 end  ----------------------------------

% 
% fprintf('总耗时：%f \n', t2) ;
   

end %idengtify_v1_3


function [RealTimeTraversing, TraversalStack] = InitRTpoint( nowTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass)
% 2.5.1 初始化遍历向量
% 初始化实时点，在开始识别或按模板边识别过程中出现没有计算到的点时使用。
% 为后续的识别提供当前模板点对应的实时点，并将生成的过程压入栈中，取栈顶过程继续后续识别。



% nowTLpoint = MarkerSetsNK(1,4) ; %第一个模板点的编号
if nowTLpoint==0 %若是识别之前的初始化
    nowTLpoint = MarkerSetsNK(1,4) ; %第一个模板点的编号
    RealTimeTraversing(idxTLpoint) = nowTLpoint ;
    RealTimeTraversing(idxTLside) = 1 ;
end

nowRT = MatchPass{nowTLpoint} ; % now_TempLate 与当前选择的模板点匹配的实时点
for j = 1:length(nowRT)
    nextRTT = RealTimeTraversing ;
    %跳过前面已经用过的点
    if ~isempty( find(nowRT(j)==nextRTT(1:nextRTT(idxTLpoint)-1), 1) ) 
        continue;
    end        
    nextRTT(nowTLpoint) = nowRT(j) ;
    TraversalStack = [TraversalStack; nextRTT] ;
end

if isempty(TraversalStack) 
    RealTimeTraversing = [] ;
    return ;
end

n = size(TraversalStack,1) ;
RealTimeTraversing = TraversalStack(n ,:) ;
TraversalStack(n ,:) = [] ;

end

function [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, MarkerSetsNK, idxTLpoint, idxTLside, ExtraStretch, MatchPass, Match_Side)
% 计算下一个模板点可能对应的实时点编号

TLside = RealTimeTraversing(idxTLside) ; %当前识别的模板边
TLMarker1 = MarkerSetsNK(TLside,3) ; %当前需要识别的模板点
nextRTT = [] ; % next_RealTimeTraversing

% 当前已识别的实时点
RTMarker2 = RealTimeTraversing(MarkerSetsNK(TLside,4)) ; 

% 若当前识别的实时点在之前没有计算到，则需要重新初始化当前实时点
if RTMarker2 == 0
    tempTLpoint = MarkerSetsNK(TLside,4) ;
    [RealTimeTraversing, TraversalStack] = InitRTpoint( tempTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass) ;
    if isempty(RealTimeTraversing), return; end
    RTMarker2 = RealTimeTraversing(MarkerSetsNK(TLside,4)) ; 
end

% 在模板边的备选边中寻找模板点Marker1可能对应的实时点
RTMarker1 = [] ;
tempRTSide = Match_Side{TLside} ;
RTMarker1 = [RTMarker1; tempRTSide(tempRTSide(:,2)==RTMarker2,3) ] ;
RTMarker1 = [RTMarker1; tempRTSide(tempRTSide(:,3)==RTMarker2,2) ] ;
RTMarker1 = unique(RTMarker1) ;

RealTimeTraversing(idxTLpoint) = TLMarker1 ;
RealTimeTraversing(idxTLside)  = TLside + 1 ;
for i = 1:length(RTMarker1)
    
    % 验证当前实时点是否已使用
    if ~isempty(find(RTMarker1(i)==RealTimeTraversing(1:TLMarker1-1),1))
        continue ;
    end
    
    % 如果当前模板点已有对应的实时点，验证之前的实时点与当前实时点是否相同
    if RealTimeTraversing(TLMarker1)>0
        if RealTimeTraversing(TLMarker1) ~= RTMarker1(i)
            continue ;
        end
    end
    
    tempRTT = RealTimeTraversing ;
    tempRTT(TLMarker1) = RTMarker1(i) ;
    nextRTT = [nextRTT; tempRTT] ;
    
end


end %cpt_nextRTT

% function plotLink(Linkage, RealTimeTraversing, xyz)
% 
% % 连线……
% 
% % Linkage               NumberOf*2    模板点之间的连接关系
% % RealTimeTraversing    1列           第i个元素表示第i个模板点对应的实时点编号
% % xyz                   n*3           重建点的3D坐标
%     
% 
% for i = 1:size(Linkage,1)
%     marker1 = RealTimeTraversing(Linkage(i,1)) ;
%     marker2 = RealTimeTraversing(Linkage(i,2)) ;
%     txyz = xyz([marker1, marker2],:) ;
%     plot3(txyz(:,1), txyz(:,2), txyz(:,3)) ;
% end
% 
% end