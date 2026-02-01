function [c] = rearrange(b,a)
c=[];
[y,index] = sort(a);
b1=sort(b);
for i=1:length(a)
    c(index(i))=b1(i);
end
end
