

function [RealTimeTraversing, flag] = identify_v2_0(MarkerSetsNK, mars, dism1m2)

% 识别
% 输入单模板的边信息、点结构体以及捕捉（采集）点的边信息，返回识别成功与否(flag)和模板点与采集点之间的识别匹配关系(RealTimeTraversing)

% v1_3  使用.mars文件中的排序方式，按模板边的顺序遍历识别，识别对象仍为模板点
% v1_4  添加模板中大于1个连通图的情况的处理
% v1_5  将识别过程中的栈 TraversalStack 设为全局变量
% v1_6  整理整个识别过程，明确变量名，去掉不需要的结构体
% v2_0  非刚体、多模板、单目标识别，模板及采集坐标信息从外部输入

% % ------------------------- 构造 MarkerSets ------------------------------
% % 用Cotex中的连线方式，连接自己匹配重建的3D坐标
% 
% % 20160906 yanjiao stand1
% MarkerSetsAll{1} = [2	1	140.696518	143.928329
% 3	1	115.64209	119.331139
% 4	2	146.049789	150.671249
% 4	3	165.742493	170.585114
% 5	1	106.658913	110.20105
% 5	2	143.60376	147.256042
% 5	3	110.913986	114.186768
% 5	4	153.373093	158.127029
% 7	6	99.973282	106.701233
% 8	6	143.855209	194.332077
% 8	7	165.964798	194.253494
% 9	6	301.714508	363.219299
% 9	7	301.026947	343.247711
% 9	8	159.819962	174.623276
% 10	8	96.866829	231.18219
% 10	9	129.771271	153.864059
% 11	9	215.815094	242.218109
% 11	10	98.9664	110.224144
% 12	11	49.058704	87.030006
% 13	11	40.200474	76.719322
% 14	6	265.746613	334.530121
% 14	7	297.130249	373.498962
% 15	6	294.679749	366.109558
% 15	7	285.185822	388.629791
% 15	14	97.378067	106.03376
% 16	14	150.934753	196.334152
% 16	15	175.273697	202.470001
% 17	14	314.366852	363.816589
% 17	15	307.643707	352.240082
% 17	16	162.387299	179.73497
% 18	16	126.203896	246.964127
% 18	17	138.161774	156.773331
% 19	17	240.821259	262.382172
% 19	18	113.630081	121.842438
% 20	19	15.791374	93.957756
% 21	19	21.806849	78.268219
% 22	1	217.066467	334.666199
% 22	2	120.241058	259.777191
% 22	3	199.743103	319.078278
% 22	4	62.843189	223.651657
% 22	6	176.130875	221.092834
% 22	7	157.651352	222.045456
% 22	14	155.262299	188.64209
% 22	15	127.55835	178.19603
% 23	22	208.193771	216.662384
% 24	23	168.756973	189.197556
% 25	24	44.258236	57.054302
% 26	6	104.588516	121.352646
% 26	14	270.213226	304.25119
% 27	23	118.697441	150.117188
% 28	26	186.879395	206.394852
% 29	14	258.508667	298.733643
% 29	28	154.973541	163.18895
% 30	28	137.147354	149.12822
% 30	29	122.473602	136.825699
% 31	6	415.31308	466.870697
% 31	28	141.215256	184.079346
% 31	30	148.016663	171.042221
% 32	7	405.252258	456.282898
% 32	25	134.248367	142.403458
% 32	31	162.584198	176.044174
% 33	14	403.692749	469.109467
% 33	29	164.796402	202.94162
% 33	30	100.575348	123.329185
% 33	31	218.460754	226.855286
% 34	15	388.955811	466.942688
% 34	25	95.1399	102.317818
% 34	27	112.552422	147.054077
% 34	32	230.771408	236.917618
% 34	33	189.051147	200.879807
% 35	31	83.290756	131.933289
% 35	32	137.469803	164.943176
% 36	31	204.910034	318.343018
% 36	32	347.556	374.109619
% 36	35	205.295288	250.515335
% 37	36	171.413071	187.916687
% 38	36	180.264053	207.873398
% 38	37	73.230843	84.237129
% 39	37	185.482452	198.357574
% 39	38	160.488373	194.606079
% 40	39	193.092102	201.109604
% 41	33	119.602417	177.742599
% 41	34	142.385147	169.778854
% 42	33	213.450043	355.936096
% 42	34	381.6875	408.223206
% 42	41	243.333786	279.772034
% 43	42	158.187302	172.685913
% 44	42	153.016373	172.195297
% 44	43	64.5224	70.541618
% 45	43	164.057465	174.112778
% 45	44	155.536118	180.993637
% 46	45	208.516693	217.349884
% ];
% 
% % % 20160717 yanjiao test12
% % MarkerSets0 = [2	1	398.127289	398.127289
% % 3	2	199.874573	199.874573
% % 4	3	400.270233	400.270233
% % 6	5	593.493408	593.493408
% % 7	6	599.131531	599.131531
% % 8	5	597.027222	597.027222
% % 8	7	590.158508	590.158508
% % ] ;
% 
% % wand500
% MarkerSetsAll{2} = [2	1	115.433189	115.433189
% 3	1	499.52655	499.52655
% 3	2	384.103271	384.103271
% ] ;
% 
% 
% % 分别输入各模板的点数、边数、阈值信息，并生成各模板（.mars）的结构体
% 
% NumberOfSide(1) = size(MarkerSetsAll{1},1) ;
% NumberOfMarker(1) = 46 ;
% useMarkerNum(1) = MarkerSets{1}(NumberOfSide(1),1) ; %模板中实际用到的Marker数量
% ExtraStretch(1) = 5 ; %误差浮动范围
% 
% NumberOfSide(2) = size(MarkerSetsAll{2},1) ;
% NumberOfMarker(2) = 3 ;
% useMarkerNum(2) = MarkerSets{2}(NumberOfSide(2),1) ; %模板中实际用到的Marker数量
% ExtraStretch(2) = 5 ; %误差浮动范围
% 
% for i = 1:length(MarkerSetsAll)
%     TemplatePointAll{i} = buildTemplatePointStruct(MarkerSetsAll{i}) ;
% end
% 
% % --------------------------- 构造MarkerSets End ----------------------------
% % ---------------------------------------------------------------------------
% 
% 
% % ---------------------------------------------------------------------------
% % -------------------------------- 整体匹配 ----------------------------------
% 
% 
% load('C:\Users\Boat\Desktop\标定匹配\input\20160906 yanjiao\nokov\huigan1_k4c3') ; %xyzn
% % load('C:\Users\Boat\Desktop\Nokov_Yanjiao_20160717\data\xyzn\xyzn_v10_5_test12') ; %xyzn
% 
% for iframe = 1:length(xyzn)
%     pause(0.015)
% % iframe = 43 ;
% % for ttttt = 1:1000
% 
% 
% flag = 0;  % 1,存在没有匹配到的模板点；2,,不存在
% xyz = xyzn{iframe} ;
% 
% % % match_v10_5_F3 之后版本生成的xyzn 需要转换
% % xyz = reshape(xyz', 3, length(xyz)/3) ;
% % xyz = xyz' ;
% 
% Numxyz = size(xyz,1) ;
% % -----------------------------------------------------------
% % 先画图……
%     markerN = size(xyz,1) ; %每帧的点数
%     fprintf('帧号：%d，点数：%d\n',iframe,markerN);
%     plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.') 
%     axis([-2300 2300 -1800 1800 -1000 2000]); % 设置坐标轴在指定的区间
%     a = [int2str(iframe),' (',int2str(markerN),')'] ;
%     title(a) ;
%     % grid on;
% % ----------------------------------------------------------
% 
% tic
% 
% % realtime中各个Marker点之间的距离
% dis = pdist(xyz,'euclidean') ; % distance，3D坐标间的欧式距离
% dis = dis' ;
% dism1m2 = zeros(length(dis),4) ; %[dis, marker1, marker2, No.#]
% dism1m2(:,1) = dis ;
% dism1m2(:,4) = (1:length(dis))' ;
% s=1 ;
% for i = 1:Numxyz-1
%     dism1m2(s:s+Numxyz-i-1,2:3) = [i*ones(Numxyz-i,1), (i+1:Numxyz)'] ;
%     s = s+Numxyz-i ;
% end
% 
% % 生成采集点之间的距离矩阵
% Nj = size(dism1m2,1) ;
% MapDis = zeros(Nj) ;
% for i = 1:Nj
%     MapDis(dism1m2(i,2), dism1m2(i,3)) = dism1m2(i,1) ;
%     MapDis(dism1m2(i,3) ,dism1m2(i,2)) = dism1m2(i,1) ;
% end
% 
% t1 = toc ;
% fprintf('计算各点距离：%f \n', t1) ;
% 
% % % 构造.mars与实时捕捉坐标的匹配结构体
% % clear MatchSets % 边的匹配关系的结构体
% %                 %  LinkageNum   模板边的编号
% %                 %  Linkage      模板边的最大最小长度以及两端点编号 
% %                 %               [MinLength, MaxLength, Marker1#, Marker2#]
% %                 %  Match        一行表示与模板边匹配的一条采集边的边长、两端点编号以及边的编号
% %                 %               [ Length, Marker1#, Marker2#, No.#]



flag = 1;  % 1,识别成功；0,识别失败
RealTimeTraversing = [] ;

% tic

NumberOfSide = mars.NumberOfSide ;
NumberOfMarker = mars.NumberOfMarker ;
useMarkerNum = mars.useMarkerNum ; %模板中实际用到的Marker数量
ExtraStretch = mars.ExtraStretch ; %误差浮动范围
TemplatePoint = mars.TemplatePoint ;

% 将采集边按模板边的长度归类，筛选模板边的备选匹配边
Match_Side = cell(NumberOfSide,1) ; % Match_TempLateSide_RealTimeSide 第i个元素存储与第i号模板边匹配的采集边

i=1; % index_MarkerSetsNK
j=1; % index_dis
Nj = size(dism1m2,1) ;
MarkerSetsNK = sortrows(MarkerSetsNK) ;
dism1m2 = sortrows(dism1m2) ;
while( i<=NumberOfSide & j<=Nj ) 
    marsNum = MarkerSetsNK(i,5) ; %模板边的编号
       
%     MatchSets(marsNum).LinkageNum = MarkerSetsNK(i,5) ;
%     MatchSets(marsNum).Linkage = MarkerSetsNK(i,1:4) ;
    
    % 寻找与当前模板边长度匹配的采集边
    if dism1m2(j,1) < MarkerSetsNK(i,1) - ExtraStretch 
        j = j+1 ;
    else
%         MatchSets(marsNum).Match = [] ;
        tempMatch = [] ;
        tj = j ;
        while (dism1m2(tj,1) < MarkerSetsNK(i,2) + ExtraStretch) 
%             MatchSets(marsNum).Match = [MatchSets(marsNum).Match; dism1m2(tj,:)] ;
            tempMatch = [tempMatch; dism1m2(tj,:)] ;
            tj = tj+1 ;
            if tj>Nj, break; end
        end
        
%         if isempty(MatchSets(marsNum).Match)
        if isempty(tempMatch)
            fprintf('%s\n',['构造边的匹配关系的结构体，模板template中第',int2str(marsNum),'号边没有匹配到~']) ;
            flag = 0 ;
            break ;
        end
        
        Match_Side{marsNum} = tempMatch ;
        i = i+1 ;
    end
end
if flag<1, return; end

t1 = toc ;
fprintf('按长度筛选备选边：%f \n', t1) ;


% 匹配一下试试看呢…… = =
Match_Point = cell(NumberOfMarker,1) ; % Match_TempLatePoint_RealTimePoint 模板点与采集点的匹配关系. 第i行存放与第i号模板点匹配的采集点编号

for i=1:NumberOfMarker
    marsNum = TemplatePoint(i).Num ;
    LinkageNum = TemplatePoint(i).LinkageNum ;
    if isempty(LinkageNum), continue; end
    
    % Step2.1 同一个模板点连接的各模板边的备选边必定含有一个公共端点
    if isempty(Match_Point{i})
%         p_re = MatchSets(LinkageNum(1)).Match(:,2:3) ; %p_re, 与模板点连接的备选边上存在的公共点
        p_re = Match_Side{LinkageNum(1)}(:,2:3) ; %p_re, 与模板点连接的备选边上存在的公共点
        p_re = unique(p_re(:)) ;
        start = 2 ;
    else
        p_re = Match_Point{i} ;
        start = 1 ;
    end
    
    for j = start:length(LinkageNum)
%         temp = MatchSets(LinkageNum(j)).Match(:,2:3) ;
        temp = Match_Side{LinkageNum(j)}(:,2:3) ;
        p_re = intersect(p_re(:),temp(:)) ;
    end
    
    if isempty(p_re) %若不存在公共点，则说明此模板点没有匹配的采集点，此帧不能整体识别
        fprintf('%s\n',['构造模板点与采集点的匹配关系，模板template中第',int2str(marsNum),'号点没有匹配到~']) ;
        flag = 0 ;
        break ;
    end
    
    Match_Point{i} = p_re ;

    % Step2.2 确定公共点为当前模板点的备选匹配点后，为另一个端点的模板点提供备选点
    LinkagePoint = TemplatePoint(i).LinkagePoint ;
    for ilp = 1:length(LinkagePoint) % index_LinkagePoint 按连接边的顺序依次比较
        
        % 寻找与该连接边另一端点匹配的采集点
%         tms = MatchSets(LinkageNum(ilp)).Match(:,2:3) ; %temp_MatchSets
        tms = Match_Side{LinkageNum(ilp)}(:,2:3) ; %temp_MatchSets
        tmp = [] ; %temp_MatchPass
        for ipre = 1:length(p_re)
            tpre = p_re(ipre) ;
            temp = tms(tms(:,1)==tpre | tms(:,2)==tpre,:) ;
            temp = temp(:) ;
            temp = temp(temp~=tpre) ;
            tmp = [tmp; temp] ; 
        end
        
        if isempty(Match_Point{LinkagePoint(ilp)})
            Match_Point{LinkagePoint(ilp)} = unique(tmp) ;
        else
            Match_Point{LinkagePoint(ilp)} = intersect(Match_Point{LinkagePoint(ilp)},tmp) ;
        end
        
        if isempty(Match_Point{LinkagePoint(ilp)})
            fprintf('%s\n',['i=',int2str(marsNum),',模板template中第',int2str(LinkagePoint(ilp)),'号点没有匹配到~']) ;
            flag = 0 ;
            break ;
        end
    end
    if flag<1, break; end
end
if flag<1, return; end


% % 生成匹配关系矩阵（同构关系矩阵）
% MapMatch = zeros(NumberOfMarker) ;
% for i = 1:size(MatchPass,1)
%     MapMatch(i,MatchPass{i}) = 1 ;
% end

t2 = toc ;
fprintf('匹配点：%f \n', t2-t1) ;
t1 = t2 ;

% % 生成real-time坐标的距离矩阵
% dis = pdist(xyz,'euclidean') ; % distance，3D坐标间的欧式距离
% Dxyz = zeros(NumberOf) ;
% for i = 1:NumberOf-1
%     Dxyz(i,i+1:NumberOf) = dis(1:NumberOf-i) ;
%     dis(1:NumberOf-i) = [] ;
% end
% Dxyz = Dxyz + Dxyz' ;


% 为模板边更新备选边
MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
RealtimeSets = [] ; %类似MarkerSets的，采集的边长与端点信息
for i = 1:NumberOfSide
%     TemplateMarker1 = MatchSets(i).Linkage(3) ;
%     TemplateMarker2 = MatchSets(i).Linkage(4) ;
    TemplateMarker1 = MarkerSetsNK(i,3) ;
    TemplateMarker2 = MarkerSetsNK(i,4) ;
    RealtimeMarker1 = Match_Point{TemplateMarker1} ;
    RealtimeMarker2 = Match_Point{TemplateMarker2} ;
    
%     tempMatch = MatchSets(i).Match ;
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
        fprintf('%s\n',['为模板边更新备选边，模板template中第',int2str(i),'号边没有匹配到~']) ;
        flag = 0 ;
        break ;
    else
%         MatchSets(i).Match = newMatch ;
        Match_Side{i} = newMatch ;
        RealtimeSets = [RealtimeSets; newMatch] ;
    end
end
if flag<1, return; end
% 去掉RealtimeSets中重复的边
 [unused1,idx,unused2] = unique(RealtimeSets(:,1)) ;
RealtimeSets = RealtimeSets(idx,:); 

t2 = toc ;
fprintf('匹配边：%f \n', t2-t1) ;
t1 = t2 ;

% % 生成采集点的邻接矩阵
% MapRealtime = zeros(useMarkerNum,size(xyz,1)) ;
% for i = 1:length(MatchSets)
%     tempMatch = MatchSets(i).Match ;
%     for j = 1:size(tempMatch,1)
%         MapRealtime(tempMatch(j,2),tempMatch(j,3)) = tempMatch(j,1) ;
%         MapRealtime(tempMatch(j,3),tempMatch(j,2)) = tempMatch(j,1)  ;
%     end
% end

% % 生成边的匹配矩阵~
% MapMatchSide = zeros(NumberOfSide, length(dis)) ;
% for i = 1:length(MatchSets)
%     LinkageNum = MatchSets(i).LinkageNum ;
%     MatchLinkageNum = MatchSets(i).Match(:,4) ;
%     MapMatchSide(LinkageNum,MatchLinkageNum) = 1 ;
% end

% =========================================================================

% 设个存档点试试看呢……
% 用模板边来匹配~？

MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
idxTLpoint = useMarkerNum+1 ; % RealTimeTraversing 中存放模板点编号的位置
idxTLside  = useMarkerNum+2 ; % RealTimeTraversing 中存放模板边编号的位置
LengthRTT  = idxTLside ;      % RealTimeTraversing 的长度
RealTimeTraversing = zeros(1,LengthRTT) ;       % 最后一个元素表示当前正在识别的模板边的编号，
                                                % 倒数第二个元素表示当前正在识别的模板点编号，
                                                % 从第一个到倒数第二个元素表示模板点对应的采集点号，
                                                % 如 RealTimeTraversing(3)=48 表示3号模板点对应的采集点号为48号
TraversalStack = [] ;

% 模板边连接的编号最小的模板点作为第一个模板点，即识别的起点，计算第一个模板点对应的采集点编号
[RealTimeTraversing, TraversalStack] = InitRTpoint( 0, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, Match_Point) ;
% nowTLpoint = MarkerSetsNK(1,4) ; %第一个模板点的编号
% RealTimeTraversing(idxTLpoint) = nowTLpoint ;
% RealTimeTraversing(idxTLside) = 1 ;
% nowRT = find( MapMatch(nowTLpoint,:)>0 )' ; % now_TempLate 与当前选择的模板点匹配的采集点
% for j = 1:length(nowRT)
%     nextRTT = RealTimeTraversing ;
%     nextRTT(nowTLpoint) = nowRT(j) ;
%     TraversalStack = [TraversalStack; nextRTT] ;
% end
% n = size(TraversalStack,1) ;
% RealTimeTraversing = TraversalStack(n ,:) ;
% TraversalStack(n ,:) = [] ;

% 按模板边的排列顺序，依次计算后续模板点对应的采集点编号
idxMarkerSets = 1 ;
while idxMarkerSets <= NumberOfSide
    
    [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, RealtimeSets, MarkerSetsNK, idxTLpoint, idxTLside, LengthRTT, ExtraStretch, Match_Point);
    
    % 入栈
    TraversalStack = [TraversalStack; nextRTT] ;
    n = size(TraversalStack,1) ;
    if n < 1
        disp('整体识别失败~') ;
        flag = 0 ;
        break ;
    end
    
    % 出栈，取栈顶的过程作为当前过程
    RealTimeTraversing = TraversalStack(n ,:) ;
    TraversalStack(n ,:) = [] ;
    idxMarkerSets = RealTimeTraversing(idxTLside) ;
end
if flag<1, return; end

RealTimeTraversing = RealTimeTraversing' ;

t2 = toc ;
fprintf('识别：%f \n', t2-t1) ;
t1 = t2 ;

fprintf('总耗时：%f \n', t2) ;
   
% hold on 
% %再连线……
% plotLink(MarkerSetsNK(:,3:4), RealTimeTraversing, xyz)
% hold off


% ttoc(ttttt) = toc ;
%     iframe
% end
% end %for iframe = 1:length(xyzn)
% mean(ttoc)
end %idengtify_v1_3


function [RealTimeTraversing, TraversalStack] = InitRTpoint( nowTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass)

% 初始化采集点，在开始识别或按模板边识别过程中出现没有计算到的点时使用。
% 为后续的识别提供当前模板点对应的采集点，并将生成的过程压入栈中，取栈顶过程继续后续识别。

% nowTLpoint = MarkerSetsNK(1,4) ; %第一个模板点的编号
if nowTLpoint==0 %若是识别之前的初始化
    nowTLpoint = MarkerSetsNK(1,4) ; %第一个模板点的编号
    RealTimeTraversing(idxTLpoint) = nowTLpoint ;
    RealTimeTraversing(idxTLside) = 1 ;
end

% nowRT = find( MapMatch(nowTLpoint,:)>0 )' ; % now_TempLate 与当前选择的模板点匹配的采集点
nowRT = MatchPass{nowTLpoint} ; % now_TempLate 与当前选择的模板点匹配的采集点
for j = 1:length(nowRT)
    nextRTT = RealTimeTraversing ;
    nextRTT(nowTLpoint) = nowRT(j) ;
    TraversalStack = [TraversalStack; nextRTT] ;
end
n = size(TraversalStack,1) ;
RealTimeTraversing = TraversalStack(n ,:) ;
TraversalStack(n ,:) = [] ;

end

function [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, RealtimeSets, MarkerSetsNK, idxTLpoint, idxTLside, LengthRTT, ExtraStretch, MatchPass)

% 计算下一个模板点可能对应的采集点编号

nowTLside = RealTimeTraversing(idxTLside) ; %当前识别的模板边
nowTLpoint = MarkerSetsNK(nowTLside,3) ; %当前识别的模板点
nextRTT = [] ;

% 当前识别的采集点
% nowRTpoint = RealTimeTraversing(RealTimeTraversing(idxTLpoint)) ; 
nowRTpoint = RealTimeTraversing(MarkerSetsNK(nowTLside,4)) ; 

% 若当前识别的采集点在之前没有计算到，则需要重新初始化当前采集点
if nowRTpoint == 0
    tempTLpoint = MarkerSetsNK(nowTLside,4) ;
    [RealTimeTraversing, TraversalStack] = InitRTpoint( tempTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass) ;
    nowRTpoint = RealTimeTraversing(MarkerSetsNK(nowTLside,4)) ; 
end

%当前采集点连接的其他采集点，以及两者间的距离
nextRTsidepoint = [RealtimeSets(RealtimeSets(:,2)==nowRTpoint,[1,3]); RealtimeSets(RealtimeSets(:,3)==nowRTpoint,[1,2])] ;

RealTimeTraversing(idxTLpoint) = nowTLpoint ;
RealTimeTraversing(idxTLside)  = nowTLside + 1 ;
for i = 1:size(nextRTsidepoint,1)
    
    % 验证当前采集点是否已使用
    if ~isempty(find(nextRTsidepoint(i,2)==RealTimeTraversing(1:nowTLpoint-1),1))
        continue ;
    end
    
    % 验证两采集点间长度是否与模板边匹配
    if nextRTsidepoint(i,1)<MarkerSetsNK(nowTLside,1)-ExtraStretch | nextRTsidepoint(i,1)>MarkerSetsNK(nowTLside,2)+ExtraStretch
        continue ;
    end
    
    % 如果当前模板点已有对应的采集点，验证之前的采集点与当前采集点是否相同
    if RealTimeTraversing(nowTLpoint)>0
        if RealTimeTraversing(nowTLpoint) ~= nextRTsidepoint(i,2)
            continue ;
        end
    end
    
    tempRTT = RealTimeTraversing ;
    tempRTT(nowTLpoint) = nextRTsidepoint(i,2) ;
    nextRTT = [nextRTT; tempRTT] ;
    
end


end %cpt_nextRTT

% function plotLink(Linkage, RealTimeTraversing, xyz)
% 
% % 连线……
% 
% % Linkage               NumberOf*2    模板点之间的连接关系
% % RealTimeTraversing    1列           第i个元素表示第i个模板点对应的采集点编号
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