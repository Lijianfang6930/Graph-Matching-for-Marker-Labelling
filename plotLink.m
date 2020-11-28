function plotLink(Linkage, RealTimeTraversing, xyz)

% 连线……

% Linkage               NumberOf*2    模板点之间的连接关系
% RealTimeTraversing    1列           第i个元素表示第i个采集点对应的模板点编号
% xyz                   n*3           重建点的3D坐标
    

for i = 1:size(Linkage,1)
    marker1 = find(RealTimeTraversing==Linkage(i,1)) ;
    marker2 = find(RealTimeTraversing==Linkage(i,2)) ;
    txyz = xyz([marker1, marker2],:) ;
    plot3(txyz(:,1), txyz(:,2), txyz(:,3)) ;
end

end