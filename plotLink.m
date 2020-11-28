function plotLink(Linkage, RealTimeTraversing, xyz)

% ���ߡ���

% Linkage               NumberOf*2    ģ���֮������ӹ�ϵ
% RealTimeTraversing    1��           ��i��Ԫ�ر�ʾ��i���ɼ����Ӧ��ģ�����
% xyz                   n*3           �ؽ����3D����
    

for i = 1:size(Linkage,1)
    marker1 = find(RealTimeTraversing==Linkage(i,1)) ;
    marker2 = find(RealTimeTraversing==Linkage(i,2)) ;
    txyz = xyz([marker1, marker2],:) ;
    plot3(txyz(:,1), txyz(:,2), txyz(:,3)) ;
end

end