

function [RealTimeTraversing, flag] = identify_v2_0(MarkerSetsNK, mars, dism1m2)

% ʶ��
% ���뵥ģ��ı���Ϣ����ṹ���Լ���׽���ɼ�����ı���Ϣ������ʶ��ɹ����(flag)��ģ�����ɼ���֮���ʶ��ƥ���ϵ(RealTimeTraversing)

% v1_3  ʹ��.mars�ļ��е�����ʽ����ģ��ߵ�˳�����ʶ��ʶ�������Ϊģ���
% v1_4  ���ģ���д���1����ͨͼ������Ĵ���
% v1_5  ��ʶ������е�ջ TraversalStack ��Ϊȫ�ֱ���
% v1_6  ��������ʶ����̣���ȷ��������ȥ������Ҫ�Ľṹ��
% v2_0  �Ǹ��塢��ģ�塢��Ŀ��ʶ��ģ�弰�ɼ�������Ϣ���ⲿ����

% % ------------------------- ���� MarkerSets ------------------------------
% % ��Cotex�е����߷�ʽ�������Լ�ƥ���ؽ���3D����
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
% % �ֱ������ģ��ĵ�������������ֵ��Ϣ�������ɸ�ģ�壨.mars���Ľṹ��
% 
% NumberOfSide(1) = size(MarkerSetsAll{1},1) ;
% NumberOfMarker(1) = 46 ;
% useMarkerNum(1) = MarkerSets{1}(NumberOfSide(1),1) ; %ģ����ʵ���õ���Marker����
% ExtraStretch(1) = 5 ; %������Χ
% 
% NumberOfSide(2) = size(MarkerSetsAll{2},1) ;
% NumberOfMarker(2) = 3 ;
% useMarkerNum(2) = MarkerSets{2}(NumberOfSide(2),1) ; %ģ����ʵ���õ���Marker����
% ExtraStretch(2) = 5 ; %������Χ
% 
% for i = 1:length(MarkerSetsAll)
%     TemplatePointAll{i} = buildTemplatePointStruct(MarkerSetsAll{i}) ;
% end
% 
% % --------------------------- ����MarkerSets End ----------------------------
% % ---------------------------------------------------------------------------
% 
% 
% % ---------------------------------------------------------------------------
% % -------------------------------- ����ƥ�� ----------------------------------
% 
% 
% load('C:\Users\Boat\Desktop\�궨ƥ��\input\20160906 yanjiao\nokov\huigan1_k4c3') ; %xyzn
% % load('C:\Users\Boat\Desktop\Nokov_Yanjiao_20160717\data\xyzn\xyzn_v10_5_test12') ; %xyzn
% 
% for iframe = 1:length(xyzn)
%     pause(0.015)
% % iframe = 43 ;
% % for ttttt = 1:1000
% 
% 
% flag = 0;  % 1,����û��ƥ�䵽��ģ��㣻2,,������
% xyz = xyzn{iframe} ;
% 
% % % match_v10_5_F3 ֮��汾���ɵ�xyzn ��Ҫת��
% % xyz = reshape(xyz', 3, length(xyz)/3) ;
% % xyz = xyz' ;
% 
% Numxyz = size(xyz,1) ;
% % -----------------------------------------------------------
% % �Ȼ�ͼ����
%     markerN = size(xyz,1) ; %ÿ֡�ĵ���
%     fprintf('֡�ţ�%d��������%d\n',iframe,markerN);
%     plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.') 
%     axis([-2300 2300 -1800 1800 -1000 2000]); % ������������ָ��������
%     a = [int2str(iframe),' (',int2str(markerN),')'] ;
%     title(a) ;
%     % grid on;
% % ----------------------------------------------------------
% 
% tic
% 
% % realtime�и���Marker��֮��ľ���
% dis = pdist(xyz,'euclidean') ; % distance��3D������ŷʽ����
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
% % ���ɲɼ���֮��ľ������
% Nj = size(dism1m2,1) ;
% MapDis = zeros(Nj) ;
% for i = 1:Nj
%     MapDis(dism1m2(i,2), dism1m2(i,3)) = dism1m2(i,1) ;
%     MapDis(dism1m2(i,3) ,dism1m2(i,2)) = dism1m2(i,1) ;
% end
% 
% t1 = toc ;
% fprintf('���������룺%f \n', t1) ;
% 
% % % ����.mars��ʵʱ��׽�����ƥ��ṹ��
% % clear MatchSets % �ߵ�ƥ���ϵ�Ľṹ��
% %                 %  LinkageNum   ģ��ߵı��
% %                 %  Linkage      ģ��ߵ������С�����Լ����˵��� 
% %                 %               [MinLength, MaxLength, Marker1#, Marker2#]
% %                 %  Match        һ�б�ʾ��ģ���ƥ���һ���ɼ��ߵı߳������˵����Լ��ߵı��
% %                 %               [ Length, Marker1#, Marker2#, No.#]



flag = 1;  % 1,ʶ��ɹ���0,ʶ��ʧ��
RealTimeTraversing = [] ;

% tic

NumberOfSide = mars.NumberOfSide ;
NumberOfMarker = mars.NumberOfMarker ;
useMarkerNum = mars.useMarkerNum ; %ģ����ʵ���õ���Marker����
ExtraStretch = mars.ExtraStretch ; %������Χ
TemplatePoint = mars.TemplatePoint ;

% ���ɼ��߰�ģ��ߵĳ��ȹ��࣬ɸѡģ��ߵı�ѡƥ���
Match_Side = cell(NumberOfSide,1) ; % Match_TempLateSide_RealTimeSide ��i��Ԫ�ش洢���i��ģ���ƥ��Ĳɼ���

i=1; % index_MarkerSetsNK
j=1; % index_dis
Nj = size(dism1m2,1) ;
MarkerSetsNK = sortrows(MarkerSetsNK) ;
dism1m2 = sortrows(dism1m2) ;
while( i<=NumberOfSide & j<=Nj ) 
    marsNum = MarkerSetsNK(i,5) ; %ģ��ߵı��
       
%     MatchSets(marsNum).LinkageNum = MarkerSetsNK(i,5) ;
%     MatchSets(marsNum).Linkage = MarkerSetsNK(i,1:4) ;
    
    % Ѱ���뵱ǰģ��߳���ƥ��Ĳɼ���
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
            fprintf('%s\n',['����ߵ�ƥ���ϵ�Ľṹ�壬ģ��template�е�',int2str(marsNum),'�ű�û��ƥ�䵽~']) ;
            flag = 0 ;
            break ;
        end
        
        Match_Side{marsNum} = tempMatch ;
        i = i+1 ;
    end
end
if flag<1, return; end

t1 = toc ;
fprintf('������ɸѡ��ѡ�ߣ�%f \n', t1) ;


% ƥ��һ�����Կ��ء��� = =
Match_Point = cell(NumberOfMarker,1) ; % Match_TempLatePoint_RealTimePoint ģ�����ɼ����ƥ���ϵ. ��i�д�����i��ģ���ƥ��Ĳɼ�����

for i=1:NumberOfMarker
    marsNum = TemplatePoint(i).Num ;
    LinkageNum = TemplatePoint(i).LinkageNum ;
    if isempty(LinkageNum), continue; end
    
    % Step2.1 ͬһ��ģ������ӵĸ�ģ��ߵı�ѡ�߱ض�����һ�������˵�
    if isempty(Match_Point{i})
%         p_re = MatchSets(LinkageNum(1)).Match(:,2:3) ; %p_re, ��ģ������ӵı�ѡ���ϴ��ڵĹ�����
        p_re = Match_Side{LinkageNum(1)}(:,2:3) ; %p_re, ��ģ������ӵı�ѡ���ϴ��ڵĹ�����
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
    
    if isempty(p_re) %�������ڹ����㣬��˵����ģ���û��ƥ��Ĳɼ��㣬��֡��������ʶ��
        fprintf('%s\n',['����ģ�����ɼ����ƥ���ϵ��ģ��template�е�',int2str(marsNum),'�ŵ�û��ƥ�䵽~']) ;
        flag = 0 ;
        break ;
    end
    
    Match_Point{i} = p_re ;

    % Step2.2 ȷ��������Ϊ��ǰģ���ı�ѡƥ����Ϊ��һ���˵��ģ����ṩ��ѡ��
    LinkagePoint = TemplatePoint(i).LinkagePoint ;
    for ilp = 1:length(LinkagePoint) % index_LinkagePoint �����ӱߵ�˳�����αȽ�
        
        % Ѱ��������ӱ���һ�˵�ƥ��Ĳɼ���
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
            fprintf('%s\n',['i=',int2str(marsNum),',ģ��template�е�',int2str(LinkagePoint(ilp)),'�ŵ�û��ƥ�䵽~']) ;
            flag = 0 ;
            break ;
        end
    end
    if flag<1, break; end
end
if flag<1, return; end


% % ����ƥ���ϵ����ͬ����ϵ����
% MapMatch = zeros(NumberOfMarker) ;
% for i = 1:size(MatchPass,1)
%     MapMatch(i,MatchPass{i}) = 1 ;
% end

t2 = toc ;
fprintf('ƥ��㣺%f \n', t2-t1) ;
t1 = t2 ;

% % ����real-time����ľ������
% dis = pdist(xyz,'euclidean') ; % distance��3D������ŷʽ����
% Dxyz = zeros(NumberOf) ;
% for i = 1:NumberOf-1
%     Dxyz(i,i+1:NumberOf) = dis(1:NumberOf-i) ;
%     dis(1:NumberOf-i) = [] ;
% end
% Dxyz = Dxyz + Dxyz' ;


% Ϊģ��߸��±�ѡ��
MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
RealtimeSets = [] ; %����MarkerSets�ģ��ɼ��ı߳���˵���Ϣ
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
        fprintf('%s\n',['Ϊģ��߸��±�ѡ�ߣ�ģ��template�е�',int2str(i),'�ű�û��ƥ�䵽~']) ;
        flag = 0 ;
        break ;
    else
%         MatchSets(i).Match = newMatch ;
        Match_Side{i} = newMatch ;
        RealtimeSets = [RealtimeSets; newMatch] ;
    end
end
if flag<1, return; end
% ȥ��RealtimeSets���ظ��ı�
 [unused1,idx,unused2] = unique(RealtimeSets(:,1)) ;
RealtimeSets = RealtimeSets(idx,:); 

t2 = toc ;
fprintf('ƥ��ߣ�%f \n', t2-t1) ;
t1 = t2 ;

% % ���ɲɼ�����ڽӾ���
% MapRealtime = zeros(useMarkerNum,size(xyz,1)) ;
% for i = 1:length(MatchSets)
%     tempMatch = MatchSets(i).Match ;
%     for j = 1:size(tempMatch,1)
%         MapRealtime(tempMatch(j,2),tempMatch(j,3)) = tempMatch(j,1) ;
%         MapRealtime(tempMatch(j,3),tempMatch(j,2)) = tempMatch(j,1)  ;
%     end
% end

% % ���ɱߵ�ƥ�����~
% MapMatchSide = zeros(NumberOfSide, length(dis)) ;
% for i = 1:length(MatchSets)
%     LinkageNum = MatchSets(i).LinkageNum ;
%     MatchLinkageNum = MatchSets(i).Match(:,4) ;
%     MapMatchSide(LinkageNum,MatchLinkageNum) = 1 ;
% end

% =========================================================================

% ����浵�����Կ��ء���
% ��ģ�����ƥ��~��

MarkerSetsNK = sortrows(MarkerSetsNK,5) ;
idxTLpoint = useMarkerNum+1 ; % RealTimeTraversing �д��ģ����ŵ�λ��
idxTLside  = useMarkerNum+2 ; % RealTimeTraversing �д��ģ��߱�ŵ�λ��
LengthRTT  = idxTLside ;      % RealTimeTraversing �ĳ���
RealTimeTraversing = zeros(1,LengthRTT) ;       % ���һ��Ԫ�ر�ʾ��ǰ����ʶ���ģ��ߵı�ţ�
                                                % �����ڶ���Ԫ�ر�ʾ��ǰ����ʶ���ģ����ţ�
                                                % �ӵ�һ���������ڶ���Ԫ�ر�ʾģ����Ӧ�Ĳɼ���ţ�
                                                % �� RealTimeTraversing(3)=48 ��ʾ3��ģ����Ӧ�Ĳɼ����Ϊ48��
TraversalStack = [] ;

% ģ������ӵı����С��ģ�����Ϊ��һ��ģ��㣬��ʶ�����㣬�����һ��ģ����Ӧ�Ĳɼ�����
[RealTimeTraversing, TraversalStack] = InitRTpoint( 0, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, Match_Point) ;
% nowTLpoint = MarkerSetsNK(1,4) ; %��һ��ģ���ı��
% RealTimeTraversing(idxTLpoint) = nowTLpoint ;
% RealTimeTraversing(idxTLside) = 1 ;
% nowRT = find( MapMatch(nowTLpoint,:)>0 )' ; % now_TempLate �뵱ǰѡ���ģ���ƥ��Ĳɼ���
% for j = 1:length(nowRT)
%     nextRTT = RealTimeTraversing ;
%     nextRTT(nowTLpoint) = nowRT(j) ;
%     TraversalStack = [TraversalStack; nextRTT] ;
% end
% n = size(TraversalStack,1) ;
% RealTimeTraversing = TraversalStack(n ,:) ;
% TraversalStack(n ,:) = [] ;

% ��ģ��ߵ�����˳�����μ������ģ����Ӧ�Ĳɼ�����
idxMarkerSets = 1 ;
while idxMarkerSets <= NumberOfSide
    
    [nextRTT, TraversalStack] = cpt_nextRTT(RealTimeTraversing, TraversalStack, RealtimeSets, MarkerSetsNK, idxTLpoint, idxTLside, LengthRTT, ExtraStretch, Match_Point);
    
    % ��ջ
    TraversalStack = [TraversalStack; nextRTT] ;
    n = size(TraversalStack,1) ;
    if n < 1
        disp('����ʶ��ʧ��~') ;
        flag = 0 ;
        break ;
    end
    
    % ��ջ��ȡջ���Ĺ�����Ϊ��ǰ����
    RealTimeTraversing = TraversalStack(n ,:) ;
    TraversalStack(n ,:) = [] ;
    idxMarkerSets = RealTimeTraversing(idxTLside) ;
end
if flag<1, return; end

RealTimeTraversing = RealTimeTraversing' ;

t2 = toc ;
fprintf('ʶ��%f \n', t2-t1) ;
t1 = t2 ;

fprintf('�ܺ�ʱ��%f \n', t2) ;
   
% hold on 
% %�����ߡ���
% plotLink(MarkerSetsNK(:,3:4), RealTimeTraversing, xyz)
% hold off


% ttoc(ttttt) = toc ;
%     iframe
% end
% end %for iframe = 1:length(xyzn)
% mean(ttoc)
end %idengtify_v1_3


function [RealTimeTraversing, TraversalStack] = InitRTpoint( nowTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass)

% ��ʼ���ɼ��㣬�ڿ�ʼʶ���ģ���ʶ������г���û�м��㵽�ĵ�ʱʹ�á�
% Ϊ������ʶ���ṩ��ǰģ����Ӧ�Ĳɼ��㣬�������ɵĹ���ѹ��ջ�У�ȡջ�����̼�������ʶ��

% nowTLpoint = MarkerSetsNK(1,4) ; %��һ��ģ���ı��
if nowTLpoint==0 %����ʶ��֮ǰ�ĳ�ʼ��
    nowTLpoint = MarkerSetsNK(1,4) ; %��һ��ģ���ı��
    RealTimeTraversing(idxTLpoint) = nowTLpoint ;
    RealTimeTraversing(idxTLside) = 1 ;
end

% nowRT = find( MapMatch(nowTLpoint,:)>0 )' ; % now_TempLate �뵱ǰѡ���ģ���ƥ��Ĳɼ���
nowRT = MatchPass{nowTLpoint} ; % now_TempLate �뵱ǰѡ���ģ���ƥ��Ĳɼ���
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

% ������һ��ģ�����ܶ�Ӧ�Ĳɼ�����

nowTLside = RealTimeTraversing(idxTLside) ; %��ǰʶ���ģ���
nowTLpoint = MarkerSetsNK(nowTLside,3) ; %��ǰʶ���ģ���
nextRTT = [] ;

% ��ǰʶ��Ĳɼ���
% nowRTpoint = RealTimeTraversing(RealTimeTraversing(idxTLpoint)) ; 
nowRTpoint = RealTimeTraversing(MarkerSetsNK(nowTLside,4)) ; 

% ����ǰʶ��Ĳɼ�����֮ǰû�м��㵽������Ҫ���³�ʼ����ǰ�ɼ���
if nowRTpoint == 0
    tempTLpoint = MarkerSetsNK(nowTLside,4) ;
    [RealTimeTraversing, TraversalStack] = InitRTpoint( tempTLpoint, RealTimeTraversing, idxTLpoint, idxTLside, TraversalStack, MarkerSetsNK, MatchPass) ;
    nowRTpoint = RealTimeTraversing(MarkerSetsNK(nowTLside,4)) ; 
end

%��ǰ�ɼ������ӵ������ɼ��㣬�Լ����߼�ľ���
nextRTsidepoint = [RealtimeSets(RealtimeSets(:,2)==nowRTpoint,[1,3]); RealtimeSets(RealtimeSets(:,3)==nowRTpoint,[1,2])] ;

RealTimeTraversing(idxTLpoint) = nowTLpoint ;
RealTimeTraversing(idxTLside)  = nowTLside + 1 ;
for i = 1:size(nextRTsidepoint,1)
    
    % ��֤��ǰ�ɼ����Ƿ���ʹ��
    if ~isempty(find(nextRTsidepoint(i,2)==RealTimeTraversing(1:nowTLpoint-1),1))
        continue ;
    end
    
    % ��֤���ɼ���䳤���Ƿ���ģ���ƥ��
    if nextRTsidepoint(i,1)<MarkerSetsNK(nowTLside,1)-ExtraStretch | nextRTsidepoint(i,1)>MarkerSetsNK(nowTLside,2)+ExtraStretch
        continue ;
    end
    
    % �����ǰģ������ж�Ӧ�Ĳɼ��㣬��֤֮ǰ�Ĳɼ����뵱ǰ�ɼ����Ƿ���ͬ
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
% % ���ߡ���
% 
% % Linkage               NumberOf*2    ģ���֮������ӹ�ϵ
% % RealTimeTraversing    1��           ��i��Ԫ�ر�ʾ��i��ģ����Ӧ�Ĳɼ�����
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