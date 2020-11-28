

function [RealTimeTraversing, flag] = identify_single_v2_6(MarkerSetsNK, mars, dism1m2)

% ʶ��
% ���뵥ģ��ı���Ϣ����ṹ���Լ���׽���ɼ�����ı���Ϣ������ʶ��ɹ����(flag)��ģ�����ʵʱ��֮���ʶ��ƥ���ϵ(RealTimeTraversing)

% identify_v1_3         ʹ��.mars�ļ��е�����ʽ����ģ��ߵ�˳�����ʶ��ʶ�������Ϊģ���
% identify_v1_4         ���ģ���д���1����ͨͼ������Ĵ���
% identify_v1_5         ��ʶ������е�ջ TraversalStack ��Ϊȫ�ֱ���
% identify_v1_6         ��������ʶ����̣���ȷ��������ȥ������Ҫ�Ľṹ��
% identify_v2_0         �Ǹ��塢��ģ�塢��Ŀ��ʶ��ģ�弰�ɼ�������Ϣ���ⲿ����
% identify_single_v2_1  ����Ϊidentify_single����ģ�嵥Ŀ��ʶ�𡣶�ģ���Ŀ����identify_wholeʵ�֡�
% identify_single_v2_2  �޸�cpt_nextRTT��Ѱ�ҵ�ǰʵʱ�����ӵ�����ʵʱ��ķ������Ż���ջ��ջ�Ĺ��̡�
% identify_single_v2_3  �ڲ��衰ȷ��������Ϊ��ǰģ���ı�ѡƥ����Ϊ��һ���˵��ģ����ṩ��ѡ�㡱��ͬʱ����ģ��ߵı�ѡ��
% identify_single_v2_4  ������ȷ��������Ϊ��ǰģ���ı�ѡƥ����Ϊ��һ���˵��ģ����ṩ��ѡ�㡱��ͬʱ����ģ��ߵı�ѡ�ߵķ���
%                       ����ģ��߷��ڸ��¶˵㱸ѡ��һ��ɾ������ RealtimeSets.
% identify_single_v2_5  ��������ಽ���У�����ģ����ģ��ߵ����ݷֿ�������ģ��ߵĲ���ŵ�ȫ��������ģ���ı�ѡ��֮��
% identify_single_v2_6  ������ģ��߹���ʱ�������������޸�ע����ʶ���ĵ��е����ƶ�Ӧ��

flag = 1;  % 1,ʶ��ɹ���0,ʶ��ʧ��
RealTimeTraversing = [] ;       % ���һ��Ԫ�ر�ʾ��ǰ����ʶ���ģ��ߵı�ţ�
                                % �����ڶ���Ԫ�ر�ʾ��ǰ����ʶ���ģ����ţ�
                                % �ӵ�һ���������ڶ���Ԫ�ر�ʾģ����Ӧ��ʵʱ��ţ�
                                % �� RealTimeTraversing(3)=48 ��ʾ3��ģ����Ӧ��ʵʱ���Ϊ48��

% tic

NumberOfSide = mars.NumberOfSide ;
NumberOfMarker = mars.NumberOfMarker ;
ExtraStretch = mars.ExtraStretch ; %������Χ
TemplatePoint = mars.TemplatePoint ; %ģ���Ľṹ��

% -----------------------  2.2 ���߳�����  --------------------------------
% ���ɼ��߰�ģ��ߵĳ��ȹ��࣬ɸѡģ��ߵı�ѡƥ���
Match_Side = cell(NumberOfSide,1) ; % Match_TempLateSide_RealTimeSide ��i��Ԫ�ش洢���i��ģ���ƥ��Ĳɼ���

i=1; % index_MarkerSetsNK
j=1; % index_dism1m2
Nj = size(dism1m2,1) ;
MarkerSetsNK = sortrows(MarkerSetsNK) ;
dism1m2 = sortrows(dism1m2) ;
while( i<=NumberOfSide  ) 
    SideNum = MarkerSetsNK(i,5) ; %ģ��ߵı��
    
    % Ѱ���뵱ǰģ��߳���ƥ��Ĳɼ���
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
%             fprintf('%s\n',['����ߵ�ƥ���ϵ�Ľṹ�壬ģ��template�е�',int2str(SideNum),'�ű�û��ƥ�䵽~']) ;
            flag = 0 ;
            break ;
        end
        
        Match_Side{SideNum} = tempMatch ;
        i = i+1 ;
    end
end
if flag<1, return; end

% t1 = toc ;
% fprintf('������ɸѡ��ѡ�ߣ�%f \n', t1) ;

% ---------------------------  2.2 end  ----------------------------------

% ------------------------  2.3 �������  ---------------------------------
% ƥ��һ�����Կ��ء��� = =
% һ��һ���ؼ���ģ���ı�ѡ��~
Match_Point = cell(NumberOfMarker,1) ; % Match_TempLatePoint_RealTimePoint ģ�����ʵʱ���ƥ���ϵ. ��i�д�����i��ģ���ƥ���ʵʱ����

for i=1:NumberOfMarker
    MarkerNum = TemplatePoint(i).Num ;
    LinkageNum = TemplatePoint(i).LinkageNum ;
    if isempty(LinkageNum), continue; end
    
    % Step2.3.1 ͬһ��ģ������ӵĸ�ģ��ߵı�ѡ�߱ض�����һ�������˵�
    if isempty(Match_Point{i})
        p_re = Match_Side{LinkageNum(1)}(:,2:3) ; %p_re, ��ģ������ӵı�ѡ���ϴ��ڵĹ�����
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
    
    if isempty(p_re) %�������ڹ����㣬��˵����ģ���û��ƥ���ʵʱ�㣬��֡��������ʶ��
%         fprintf('%s\n',['����ģ�����ʵʱ���ƥ���ϵ��ģ��template�е�',int2str(MarkerNum),'�ŵ�û��ƥ�䵽~']) ;
        flag = 0 ;
        break ;
    end
    
    Match_Point{i} = p_re ;

    % Step2.3.2 ȷ��������Ϊ��ǰģ���ı�ѡƥ����Ϊ��һ���˵��ģ����ṩ��ѡ��
    LinkagePoint = TemplatePoint(i).LinkagePoint ;
    for ilp = 1:length(LinkagePoint) % index_LinkagePoint �����ӱߵ�˳�����αȽ�
        
        % Ѱ��������ӱ���һ�˵�ƥ���ʵʱ��Ͷ�Ӧ�Ĳɼ���
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
        
        % ���µ�ͱߵ�ƥ����Ϣ
        % ����������ĵ㡢����ƥ����еĵ㡢��ƥ����Ϣ�Ƚϣ��������߶����ڵ���Ŀ
        if isempty(Match_Point{LinkagePoint(ilp)})
            Match_Point{LinkagePoint(ilp)} = nowMP ;
        else
            newMP = intersect(Match_Point{LinkagePoint(ilp)},nowMP) ;   %new_Match_Point
            if isempty(newMP)
%                 fprintf('%s\n',['i=',int2str(MarkerNum),',ģ��template�е�',int2str(LinkagePoint(ilp)),'�ŵ�û��ƥ�䵽~']) ;
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

% --------------------  2.4 ����ģ��߶�Ӧ�ı�ѡ��  ----------------------------

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
%         fprintf('%s\n',['Ϊģ��߸��±�ѡ�ߣ�ģ��template�е�',int2str(i),'�ű�û��ƥ�䵽~']) ;
        flag = 0 ;
        break ;
    else
        Match_Side{i} = newMatch ;
    end
end
if flag<1, return; end

% ---------------------------  2.4 end  ----------------------------------

% --------------------  2.5 ��ģ���˳�����ʶ��  --------------------------

% ����浵�����Կ��ء���
% ��ģ�����ƥ��~��

MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
idxTLpoint = NumberOfMarker+1 ; % RealTimeTraversing �д��ģ����ŵ�λ��
idxTLside  = NumberOfMarker+2 ; % RealTimeTraversing �д��ģ��߱�ŵ�λ��
LengthRTT  = idxTLside ;      % RealTimeTraversing �ĳ���
RealTimeTraversing = zeros(1,LengthRTT) ;       % ���һ��Ԫ�ر�ʾ��ǰ����ʶ���ģ��ߵı�ţ�
                                                % �����ڶ���Ԫ�ر�ʾ��ǰ����ʶ���ģ����ţ�
                                                % �ӵ�һ���������ڶ���Ԫ�ر�ʾģ����Ӧ��ʵʱ��ţ�
                                                % �� RealTimeTraversing(3)=48 ��ʾ3��ģ����Ӧ��ʵʱ���Ϊ48��
TraversalStack = [] ;

% 2.5.1 ��ʼ����������
% ģ������ӵı����С��ģ�����Ϊ��һ��ģ��㣬��ʶ�����㣬�����һ��ģ����Ӧ��ʵʱ����
[RealTimeTraversing, TraversalStack] = InitRTpoint( 0, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, Match_Point) ;

% 2.5.2 ������һ��ģ����ʶ���
% ��ģ��ߵ�����˳�����μ������ģ����Ӧ��ʵʱ����
idxMarkerSets = 1 ;
while idxMarkerSets <= NumberOfSide
    
    % ������һ�����ܳ��ֵ�RealTimeTraversing
    [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, MarkerSetsNK, idxTLpoint, idxTLside, ExtraStretch, Match_Point, Match_Side);
    
    nRTT = size(nextRTT,1) ; %nextRTT�к��еĹ�����
    if nRTT == 1
        RealTimeTraversing = nextRTT ;
    elseif nRTT < 1
        n = size(TraversalStack,1) ;
        if n < 1
%             disp('����ʶ��ʧ��~') ;
            flag = 0 ;
            break ;
        end
        % ��ջ��ȡջ���Ĺ�����Ϊ��ǰ����
        RealTimeTraversing = TraversalStack(n ,:) ;
        TraversalStack(n ,:) = [] ;     
    else
        RealTimeTraversing = nextRTT(nRTT,:) ;
        % ��ջ
        TraversalStack = [TraversalStack; nextRTT(1:nRTT-1,:)] ;
    end
    idxMarkerSets = RealTimeTraversing(idxTLside) ;
end
if flag<1, return; end

RealTimeTraversing = RealTimeTraversing' ;

% t2 = toc ;
% fprintf('ʶ��%f \n', t2-t1) ;
% t1 = t2 ;

% ---------------------------  2.5 end  ----------------------------------

% 
% fprintf('�ܺ�ʱ��%f \n', t2) ;
   

end %idengtify_v1_3


function [RealTimeTraversing, TraversalStack] = InitRTpoint( nowTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass)
% 2.5.1 ��ʼ����������
% ��ʼ��ʵʱ�㣬�ڿ�ʼʶ���ģ���ʶ������г���û�м��㵽�ĵ�ʱʹ�á�
% Ϊ������ʶ���ṩ��ǰģ����Ӧ��ʵʱ�㣬�������ɵĹ���ѹ��ջ�У�ȡջ�����̼�������ʶ��



% nowTLpoint = MarkerSetsNK(1,4) ; %��һ��ģ���ı��
if nowTLpoint==0 %����ʶ��֮ǰ�ĳ�ʼ��
    nowTLpoint = MarkerSetsNK(1,4) ; %��һ��ģ���ı��
    RealTimeTraversing(idxTLpoint) = nowTLpoint ;
    RealTimeTraversing(idxTLside) = 1 ;
end

nowRT = MatchPass{nowTLpoint} ; % now_TempLate �뵱ǰѡ���ģ���ƥ���ʵʱ��
for j = 1:length(nowRT)
    nextRTT = RealTimeTraversing ;
    %����ǰ���Ѿ��ù��ĵ�
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
% ������һ��ģ�����ܶ�Ӧ��ʵʱ����

TLside = RealTimeTraversing(idxTLside) ; %��ǰʶ���ģ���
TLMarker1 = MarkerSetsNK(TLside,3) ; %��ǰ��Ҫʶ���ģ���
nextRTT = [] ; % next_RealTimeTraversing

% ��ǰ��ʶ���ʵʱ��
RTMarker2 = RealTimeTraversing(MarkerSetsNK(TLside,4)) ; 

% ����ǰʶ���ʵʱ����֮ǰû�м��㵽������Ҫ���³�ʼ����ǰʵʱ��
if RTMarker2 == 0
    tempTLpoint = MarkerSetsNK(TLside,4) ;
    [RealTimeTraversing, TraversalStack] = InitRTpoint( tempTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass) ;
    if isempty(RealTimeTraversing), return; end
    RTMarker2 = RealTimeTraversing(MarkerSetsNK(TLside,4)) ; 
end

% ��ģ��ߵı�ѡ����Ѱ��ģ���Marker1���ܶ�Ӧ��ʵʱ��
RTMarker1 = [] ;
tempRTSide = Match_Side{TLside} ;
RTMarker1 = [RTMarker1; tempRTSide(tempRTSide(:,2)==RTMarker2,3) ] ;
RTMarker1 = [RTMarker1; tempRTSide(tempRTSide(:,3)==RTMarker2,2) ] ;
RTMarker1 = unique(RTMarker1) ;

RealTimeTraversing(idxTLpoint) = TLMarker1 ;
RealTimeTraversing(idxTLside)  = TLside + 1 ;
for i = 1:length(RTMarker1)
    
    % ��֤��ǰʵʱ���Ƿ���ʹ��
    if ~isempty(find(RTMarker1(i)==RealTimeTraversing(1:TLMarker1-1),1))
        continue ;
    end
    
    % �����ǰģ������ж�Ӧ��ʵʱ�㣬��֤֮ǰ��ʵʱ���뵱ǰʵʱ���Ƿ���ͬ
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
% % ���ߡ���
% 
% % Linkage               NumberOf*2    ģ���֮������ӹ�ϵ
% % RealTimeTraversing    1��           ��i��Ԫ�ر�ʾ��i��ģ����Ӧ��ʵʱ����
% % xyz                   n*3           �ؽ����3D����
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