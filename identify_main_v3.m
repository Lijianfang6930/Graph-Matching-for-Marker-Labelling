

function time = identify_main_v3(xyzn)

% ʶ��

% identify_v1_3     ʹ��.mars�ļ��е�����ʽ����ģ��ߵ�˳�����ʶ��ʶ�������Ϊģ���
% identify_v1_4     ���ģ���д���1����ͨͼ������Ĵ���
% identify_v1_5     ��ʶ������е�ջ TraversalStack ��Ϊȫ�ֱ���
% identify_v1_6     ��������ʶ����̣���ȷ��������ȥ������Ҫ�Ľṹ��
% identify_v1_7     �Ǹ��塢��ģ�塢��Ŀ��ʶ��
% identify_main_v1  ����Ϊidentify_main���ֶ�����ɼ��ĸ�֡�����Լ�ģ����Ϣ��
%                   ��������ʶ����identify_whole����ģ�塢��Ŀ�����ʶ��
% identify_main_v2  ��Ӹ��ٵļ�����
% identify_main_v3  identify_whole�н�һ��ģ����Զ�Ӧ�����ͬĿ��ķ�������Ϊһ��ģ��ֻ��Ӧһ��Ŀ��
%                   һ��ģ���Ӧ�����ͬĿ��ķ�������������

% ------------------------- ���� MarkerSets ------------------------------
% ��Cotex�е����߷�ʽ�������Լ�ƥ���ؽ���3D����

% % 20160906 yanjiao stand1
MarkerSetsAll{1} = [2	1	140.696518	143.928329
3	1	115.64209	119.331139
4	2	146.049789	150.671249
4	3	165.742493	170.585114
5	1	106.658913	110.20105
5	2	143.60376	147.256042
5	3	110.913986	114.186768
5	4	153.373093	158.127029
7	6	99.973282	106.701233
8	6	143.855209	194.332077
8	7	165.964798	194.253494
9	6	301.714508	363.219299
9	7	301.026947	343.247711
9	8	159.819962	174.623276
10	8	96.866829	231.18219
10	9	129.771271	153.864059
11	9	215.815094	242.218109
11	10	98.9664	110.224144
12	11	49.058704	87.030006
13	11	40.200474	76.719322
14	6	265.746613	334.530121
14	7	297.130249	373.498962
15	6	294.679749	366.109558
15	7	285.185822	388.629791
15	14	97.378067	106.03376
16	14	150.934753	196.334152
16	15	175.273697	202.470001
17	14	314.366852	363.816589
17	15	307.643707	352.240082
17	16	162.387299	179.73497
18	16	126.203896	246.964127
18	17	138.161774	156.773331
19	17	240.821259	262.382172
19	18	113.630081	121.842438
20	19	15.791374	93.957756
21	19	21.806849	78.268219
22	1	217.066467	334.666199
22	2	120.241058	259.777191
22	3	199.743103	319.078278
22	4	62.843189	223.651657
22	6	176.130875	221.092834
22	7	157.651352	222.045456
22	14	155.262299	188.64209
22	15	127.55835	178.19603
23	22	208.193771	216.662384
24	23	168.756973	189.197556
25	24	44.258236	57.054302
26	6	104.588516	121.352646
26	14	270.213226	304.25119
27	23	118.697441	150.117188
28	26	186.879395	206.394852
29	14	258.508667	298.733643
29	28	154.973541	163.18895
30	28	137.147354	149.12822
30	29	122.473602	136.825699
31	6	415.31308	466.870697
31	28	141.215256	184.079346
31	30	148.016663	171.042221
32	7	405.252258	456.282898
32	25	134.248367	142.403458
32	31	162.584198	176.044174
33	14	403.692749	469.109467
33	29	164.796402	202.94162
33	30	100.575348	123.329185
33	31	218.460754	226.855286
34	15	388.955811	466.942688
34	25	95.1399	102.317818
34	27	112.552422	147.054077
34	32	230.771408	236.917618
34	33	189.051147	200.879807
35	31	83.290756	131.933289
35	32	137.469803	164.943176
36	31	204.910034	318.343018
36	32	347.556	374.109619
36	35	205.295288	250.515335
37	36	171.413071	187.916687
38	36	180.264053	207.873398
38	37	73.230843	84.237129
39	37	185.482452	198.357574
39	38	160.488373	194.606079
40	39	193.092102	201.109604
41	33	119.602417	177.742599
41	34	142.385147	169.778854
42	33	213.450043	355.936096
42	34	381.6875	408.223206
42	41	243.333786	279.772034
43	42	158.187302	172.685913
44	42	153.016373	172.195297
44	43	64.5224	70.541618
45	43	164.057465	174.112778
45	44	155.536118	180.993637
46	45	208.516693	217.349884
];

% % 20160717 yanjiao test12
% MarkerSets0 = [2	1	398.127289	398.127289
% 3	2	199.874573	199.874573
% 4	3	400.270233	400.270233
% 6	5	593.493408	593.493408
% 7	6	599.131531	599.131531
% 8	5	597.027222	597.027222
% 8	7	590.158508	590.158508
% ] ;

% % Lframe
% MarkerSetsAll{1} = [2	1	401.002991	401.640076
% 3	1	447.144745	448.059235
% 3	2	198.007111	199.51741
% 4	1	721.802246	722.687927
% 4	2	598.112549	599.316162
% 4	3	399.597656	400.236877
% ] ;

% % wand500
% MarkerSetsAll{2} = [2	1	115.433189	115.433189
% 3	1	499.52655	499.52655
% 3	2	384.103271	384.103271
% ] ;


% Board4_1
MarkerSetsAll{5} = [3	1	127.539803	127.539803
3	2	94.090744	94.090744
4	1	95.530128	95.530128
4	2	129.16069	129.16069
4	3	83.988762	83.988762
];

% Board4_2
MarkerSetsAll{6} = [3	1	109.13018	109.13018
3	2	94.717377	94.717377
4	2	128.811508	128.811508
4	3	62.340546	62.340546
];

% Board4_3
MarkerSetsAll{7} = [2	1	48.988018	48.988018
3	1	88.516998	88.516998
3	2	93.696907	93.696907
4	1	60.011986	60.011986
4	3	63.965576	63.965576
];

% Board4_4
MarkerSetsAll{8} = [2	1	63.087414	63.087414
3	2	60.096985	60.096985
4	1	84.052063	84.052063
4	2	66.458603	66.458603
4	3	92.066238	92.066238
];

% Board5_1
MarkerSetsAll{1} = [2	1	38.769585	38.769585
3	1	63.979641	63.979641
3	2	47.45274	47.45274
4	1	106.194672	106.194672
4	3	50.234104	50.234104
5	2	109.55619	109.55619
5	3	67.878242	67.878242
5	4	46.6679	46.6679
];

% Board5_2
MarkerSetsAll{2} = [3	1	104.380554	104.380554
4	1	100.463821	100.463821
4	2	87.065544	87.065544
5	2	94.153313	94.153313
5	3	103.449524	103.449524
];

% Board5_3
MarkerSetsAll{3} = [2	1	59.228962	59.228962
3	2	37.060242	37.060242
4	1	83.460823	83.460823
4	2	103.123581	103.123581
5	1	102.847069	102.847069
5	2	84.505707	84.505707
5	3	91.425636	91.425636
5	4	59.990185	59.990185
];

% Board5_4
MarkerSetsAll{4} = [2	1	86.556702	86.556702
3	1	62.043701	62.043701
4	3	66.286629	66.286629
5	3	66.469116	66.469116
5	4	86.897446	86.897446
];


% �ֱ������ģ��ĵ�������������ֵ��Ϣ�������ɸ�ģ�壨.mars���Ľṹ��
% mars
% NumberOfSide      ģ���бߵ�������
% NumberOfMarker    ģ���е���ܸ���
% ExtraStretch      ��ƥ��ʱ����ֵ

for i=1:4
    mars(i).NumberOfSide = size(MarkerSetsAll{i},1) ;
    mars(i).NumberOfMarker = 5 ;
    mars(i).ExtraStretch = 3 ; 
    mars(i+4).NumberOfSide = size(MarkerSetsAll{i+4},1) ;
    mars(i+4).NumberOfMarker = 4 ;
    mars(i+4).ExtraStretch = 3 ; 
end
% 
% mars(1).NumberOfSide = size(MarkerSetsAll{1},1) ;
% mars(1).NumberOfMarker = 46 ;
% mars(1).ExtraStretch = 20 ; %������Χ
% 
% % mars(2).NumberOfSide = size(MarkerSetsAll{2},1) ;
% % mars(2).NumberOfMarker = 3 ;
% % mars(2).ExtraStretch = 8 ; %������Χ

for i = 1:length(MarkerSetsAll)
    mars(i).TemplatePoint = buildTemplatePointStruct(MarkerSetsAll{i}, mars(i).NumberOfMarker) ;
    tempMarkerSets = MarkerSetsAll{i} ;
    tempMarkerSets=[tempMarkerSets(:,3:4), tempMarkerSets(:,1:2), (1:size(tempMarkerSets,1))'] ;%[ MinLength, MaxLength, Marker1#, Marker2#, No.# ]
    MarkerSetsAll{i} = tempMarkerSets ;
end

% --------------------------- ����MarkerSets End ----------------------------
% ---------------------------------------------------------------------------


% ---------------------------------------------------------------------------
% -------------------------------- ����ƥ�� ----------------------------------

% ����ɢ������
% load('C:\Users\Boat\Desktop\�궨ƥ��\input\20160906 yanjiao\nokov\xiadun1_k4c3') ; %xyzn
% load('C:\Users\Boat\Desktop\Nokov_Yanjiao_20160717\data\xyzn\xyzn_v10_5_test12') ; %xyzn
tempxyz = importdata(['input\Cotex\8Board1.txt']) ;
for i = 1:length(tempxyz)
    xyzn{i} = tempxyz(i,:) ;
end


MarsUsed = zeros(1,length(mars)) ; % һ�У���i��ֵ��ʾ��i��ģ��ĵ�ǰʶ�������1,ʶ����Ŀ�ꣻ0,û��ʶ��
Traversing_mars = cell(1,length(mars)) ; % ��i��cell�洢��i��ģ��ʶ�������Ŀ�ꣻ
                                         % ÿ��cell�У�һ�б�ʾһ��Ŀ���ʶ��㣬��i�б�ʾ��i��ģ����Ӧ�Ĳɼ���(xyz0)�ı��
xyz_indexLast = [] ; %��һ֡�Ĳɼ������꣬������ʹ��
         

time = nan(1,length(xyzn)) ;

for iframe = 1:1:length(xyzn)
      %iframe = 301 ; %debug
    pause(0.015)
% iframe = 43 ;
% for ttttt = 1:1000
tic

xyz = xyzn{iframe} ;
Numxyz = size(xyz,1) ; %�ɼ��������������match_v10_5_F3֮ǰ�汾���ɵ�xyzn

% match_v10_5_F3 ֮��汾���ɵ�xyzn ��Ҫת��
Numxyz = length(xyz)/3 ;
xyz = reshape(xyz', 3, Numxyz) ;
xyz = xyz' ;

xyz_index = [xyz , (1:Numxyz)'] ; %���ɼ�����
xyz_index0 = xyz_index ; %�ݴ��֡�Ĳɼ��㼰���

% -----------------------------------------------------------
% ��ɢ�㻭ͼ����
    markerN = size(xyz,1) ; %ÿ֡�ĵ���
%     fprintf('֡�ţ�%d��������%d\n',iframe,markerN);
    plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.') 
%     axis([-2300 2300 -1800 1800 -1000 2000]); % ������������ָ��������
%     axis([-1800 1800 -1300 1300 -500 1500]); % ������������ָ��������
    axis([-1500 1500 -1000 1000 0 1400]); % ������������ָ��������
    a = [int2str(iframe),' (',int2str(markerN),')'] ;
    title(a) ;
    % grid on;
% ----------------------------------------------------------

% if iframe==67
%     iframe %debug
% end

% ����
% [Traversing_mars, MarsUsed, xyz_index] = tracking(xyz_index, xyz_indexLast, MarkerSetsAll, mars, Traversing_mars, MarsUsed) ; 


% ����ʶ��
MarsUsed = zeros(1,length(mars)) ;  %debug
Traversing_mars = cell(1,length(mars)) ; %debug
[Traversing_mars, MarsUsed] = identify_whole(xyz_index, MarkerSetsAll, mars, Traversing_mars, MarsUsed) ;


% ���ߡ���
plotLink(MarkerSetsAll, Traversing_mars, xyz, mars) ;

%Ϊ��һ֡�ĸ����ṩ����
xyz_indexLast = xyz_index0 ; 

time(1,iframe) = toc ;
time(2,iframe) = ( sum(MarsUsed) == length(MarsUsed) ) ; %����ʶ��ȫ��ģ��

end %for iframe = 1:length(xyzn)

end %idengtify_main


function TemplatePoint = buildTemplatePointStruct(MarkerSets, NumberOfMarker)

% ����.mars��ģ��ߵ���Ϣ������ģ���Ľṹ��
% MarkerSets    4��              [Marker1#, Marker2#, MinLength, MaxLength]
% TemplatePoint ģ���Ľṹ��     Num           No.#     ģ���ı��
%                                 LinkageNum    һ��     ��ǰ�����ڵ�ģ��ߵı��
%                                 LinkagePoint  һ��     �뵱ǰ��������ģ���ı��
%                                 Dis           ����     ��ǰ�����ڵ�ģ��ߵ���С����󳤶�

MarkerSets=[MarkerSets(:,3:4), MarkerSets(:,1:2), (1:size(MarkerSets,1))'] ;%[ MinLength, MaxLength, Marker1#, Marker2#, No.# ]

TemplatePoint = [] ;   
for i = 1:NumberOfMarker
    TemplatePoint(i).Num = i ;
    idx = MarkerSets(:,3)==i ;
    TemplatePoint(i).LinkageNum = MarkerSets(idx,5) ;
    TemplatePoint(i).LinkagePoint = MarkerSets(idx,4) ;
    TemplatePoint(i).Dis  = MarkerSets(idx,1:2) ;
    idx = MarkerSets(:,4)==i ;
    TemplatePoint(i).LinkageNum = [TemplatePoint(i).LinkageNum; MarkerSets(idx,5)] ;
    TemplatePoint(i).LinkagePoint = [TemplatePoint(i).LinkagePoint; MarkerSets(idx,3)] ;
    TemplatePoint(i).Dis  = [TemplatePoint(i).Dis; MarkerSets(idx,1:2)] ;
end

end %buildTemplatePointStruct


function [Traversing_mars, MarsUsed] = identify_whole(xyz_index, MarkerSetsAll, mars, Traversing_mars, MarsUsed) 

% ����ʶ��
% û������ʶ�𵽵�ģ���������ʶ����ʶ�𵽵�ģ�岻������


% Numxyz = size(xyz_index,1) ;
% xyz0 = [xyz , (1:Numxyz)'] ; %���ɼ�����

% realtime�и���Marker��֮��ľ���
xyz = xyz_index(:,1:3) ;
dis = pdist(xyz,'euclidean') ; % distance��3D������ŷʽ����
dis = dis' ;
dism1m2 = zeros(length(dis),4) ; %[dis, marker1, marker2, No.#]
dism1m2(:,1) = dis ;
dism1m2(:,4) = (1:length(dis))' ;
s=1 ;
Numxyz = size(xyz,1) ;
for i = 1:Numxyz-1
    dism1m2(s:s+Numxyz-i-1,2:3) = [i*ones(Numxyz-i,1), (i+1:Numxyz)'] ;
    s = s+Numxyz-i ;
end

% һ��ģ��һ��ģ���ʶ�𡭡�
missmars = find(MarsUsed<1) ;
RealTimeTraversing = [] ;       % ���һ��Ԫ�ر�ʾ��ǰ����ʶ���ģ��ߵı�ţ�
                                % �����ڶ���Ԫ�ر�ʾ��ǰ����ʶ���ģ����ţ�
                                % �ӵ�һ���������ڶ���Ԫ�ر�ʾģ����Ӧ��ʵʱ��ţ�
                                % �� RealTimeTraversing(3)=48 ��ʾ3��ģ����Ӧ��ʵʱ���Ϊ48��
flag = 0;
for imars = missmars
    
%     fprintf('ģ��%d\n ',imars) ;
%   tic  
        %��ǰһ��ģ��ʶ��ɹ�������ʶ��ĵ�ȥ�������������ļ���
        
        xyz_index(RealTimeTraversing(1:length(RealTimeTraversing)-2), :) = [] ; %ȥ��ǰ��ģ���Ѿ�ʶ������ĵ�
        
        %��ʣ��ɼ����������ģ��������٣��������ǰģ���ʶ��
        if size(xyz_index,1) < mars(imars).NumberOfMarker, break; end
        
        % realtime�и���Marker��֮��ľ���
        xyz = xyz_index(:,1:3) ;
        dis = pdist(xyz,'euclidean') ; % distance��3D������ŷʽ����
        dis = dis' ;
        dism1m2 = zeros(length(dis),4) ; %[dis, marker1, marker2, No.#]
        dism1m2(:,1) = dis ;
        dism1m2(:,4) = (1:length(dis))' ;
        s=1 ;
        Numxyz = size(xyz,1) ;
        for i = 1:Numxyz-1
            dism1m2(s:s+Numxyz-i-1,2:3) = [i*ones(Numxyz-i,1), (i+1:Numxyz)'] ;
            s = s+Numxyz-i ;
        end
        
% t1 = toc ;
% fprintf('���������룺%f \n', t1) ;
        
        % ��ģ�嵥Ŀ��ʶ��
        [RealTimeTraversing, flag] = identify_single_v2_6(MarkerSetsAll{imars}, mars(imars), dism1m2) ;
        
        % ��ʶ��ɹ��������ģ���Ӧ��ʶ�����, flag 1,ʶ��ɹ���0,ʶ��ʧ��
        if flag>0
%             fprintf('ʶ��ɹ�\n') ;
            tempTraversing = RealTimeTraversing ;
            for iTT = 1 : mars(imars).NumberOfMarker
                tempTraversing(iTT) = xyz_index( tempTraversing(iTT), 4) ;
            end
            
            Traversing_mars{imars} =  tempTraversing(1:mars(imars).NumberOfMarker)' ;
            MarsUsed(imars) = 1 ;
            
%             hold on
%             %�����ߡ���
%             plotLink_single(MarkerSetsAll{imars}(:,3:4), RealTimeTraversing, xyz)
%             hold off
        else
            RealTimeTraversing = [] ;
        end
        
end %for imars = 1:length(mars)

end %Traversing_mars


function [Traversing_mars, MarsUsed, xyz_index] = tracking(xyz_index, xyz_indexLast, MarkerSetsAll, mars, Traversing_mars, MarsUsed) 

% ����

Threshold = 50 ;

% һ��ģ��һ��ģ��ظ��١���
for imars = find(MarsUsed>0)
    
    % һ��Ŀ��һ��Ŀ��ظ��١���
    iopject = 1 ;
    while iopject <= size(Traversing_mars{imars},1)
        Traversing = Traversing_mars{imars}(iopject,:) ;
        
        % һ����һ����ظ��١���
        for i = 1:length(Traversing)
            if Traversing(i)==0, continue; end
            
            % Ԥ���i�����ڵ�ǰ֡�е�λ�÷�Χ
            XyzPointLast = xyz_indexLast(Traversing(i),1:3) ;
            Distance = bsxfun(@minus, xyz_index(:,1:3), XyzPointLast ) ; % distance = data - repmat(data(idx,:),size(data,1),1) ;
            Distance = sum(abs(Distance),2) ;
            idx = find(Distance < Threshold) ; %���ٵ��ĵ���±�
            
            if length(idx) < 1 %û�и��ٵ���
%                 fprintf('ģ��ţ�%d��Ŀ��ţ�%d��ģ���ţ�%d�����ٵ��ĵ�����%d\n',imars,iopject,i,length(idx)) ;
                Traversing(i) = 0 ;
            elseif length(idx) > 1 %���ٵ�����1����
%                 fprintf('ģ��ţ�%d��Ŀ��ţ�%d��ģ���ţ�%d�����ٵ��ĵ�����%d\n',imars,iopject,i,length(idx)) ;
                Traversing(i) = xyz_index(idx(1),4) ; %�����ȡһ�����԰ɡ��� = =
                xyz_index(idx,:) = [] ; %ɾ���Ѿ����ٵ��ĵ㣬������֮��ĸ��١�ʶ�����                
            else %���ٵ�Ψһ��
                Traversing(i) = xyz_index(idx,4) ;
                xyz_index(idx,:) = [] ; %ɾ���Ѿ����ٵ��ĵ㣬������֮��ĸ��١�ʶ�����
            end
            
        end
        
        % ����Ŀ������˰��������ϣ�����Ϊ��Ŀ�����ʧ�ܣ���ʶ��ģ����ɾ����Ŀ��
        if sum(Traversing<1) > length(Traversing)/2
            Traversing_mars{imars}(iopject,:) = [] ;
        else
            Traversing_mars{imars}(iopject,:) = Traversing ;
            iopject = iopject + 1;
        end
        
    end % for iopject <= size(Traversing_mars{imars},1)
    
    % ��һ��ģ��δ���ٵ��κ�Ŀ�꣬��MarsUsed�н���Ŀ����Ϊδʶ��״̬
    if isempty(Traversing_mars{imars})
        MarsUsed(imars) = 0 ;
    end
    
end % for imars = find(MarsUsed>0)

end



% function plotLink_single(Linkage, RealTimeTraversing, xyz)
% 
% % ���ߡ���
% % ���ڵ�ģ�嵥Ŀ��ʶ��
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

function plotLink(MarkerSetsAll, Traversing_mars, xyz, mars) 

% ���ߡ���
% ����һ֡�Ķ�ģ���Ŀ��ʶ��

% Linkage               NumberOf*2    ģ���֮������ӹ�ϵ
% RealTimeTraversing    1��           ��i��Ԫ�ر�ʾ��i��ģ����Ӧ�Ĳɼ�����
% xyz                   n*3           �ؽ����3D����

hold on

% һ��Ŀ��һ��Ŀ���������
for imars = 1:length(Traversing_mars) 
    Linkage = MarkerSetsAll{imars}(:,3:4) ;
    for iopject = 1:size(Traversing_mars{imars},1)
        RealTimeTraversing = Traversing_mars{imars}(iopject,:) ;
        
        for i = 1:mars(imars).NumberOfSide
            marker1 = RealTimeTraversing(Linkage(i,1)) ;
            marker2 = RealTimeTraversing(Linkage(i,2)) ;
            if marker1==0 | marker2==0, continue; end
            txyz = xyz([marker1, marker2],:) ;
            plot3(txyz(:,1), txyz(:,2), txyz(:,3)) ;
        end
    end
end

hold off

end

