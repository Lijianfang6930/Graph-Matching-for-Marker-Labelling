function b = Repeated(p)
%��ȡp���ظ���Ԫ��
% input
% p      ���������
% output
% b      1��  p���ظ���Ԫ��

p=p(:) ;
b = sort(p);
db = diff(b);
b = b(db==0) ; %�ظ����ֵ�Ԫ��
if isempty(b)
    return ;
end
db = diff(b);
d = db ~= 0;
d(numel(b)) = true; % Final element is always a member of unique list.
b = b(d);

end