    
function RealTimeTraversing = SelectSide(RealTimeTraversing, nowpoint, RealtimeSets, TraversalStack, StackLength)

%     nowpoint = now(length(now)) ;
nowside = RealtimeSets( RealtimeSets(:,2)==nowpoint | RealtimeSets(:,3)==nowpoint , 4 ) ;

if length(nowside)>1 
    for j = 1:length(nowside)-1
        TraversalStack{StackLength} = [RealTimeTraversing, nowside(j)] ;
        StackLength = StackLength + 1 ;
    end
    nowside = nowside(j+1) ;
    RealTimeTraversing = [RealTimeTraversing, nowside] ;
    visitedSide(nowside) = 1 ;
elseif length(nowside)<1 
    
    
    
else
    RealTimeTraversing = [RealTimeTraversing, nowside] ;
    visitedSide(nowside) = 1 ;
end

