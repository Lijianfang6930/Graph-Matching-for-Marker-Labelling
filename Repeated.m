function b = Repeated(p)
%提取p中重复的元素
% input
% p      矩阵或向量
% output
% b      1行  p中重复的元素

p=p(:) ;
b = sort(p);
db = diff(b);
b = b(db==0) ; %重复出现的元素
if isempty(b)
    return ;
end
db = diff(b);
d = db ~= 0;
d(numel(b)) = true; % Final element is always a member of unique list.
b = b(d);

end