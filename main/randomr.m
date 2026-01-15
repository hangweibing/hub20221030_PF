function [outindex] = randomr(w)
N=length(w);
parfor i=1:N
    outindex(i)=find(rand<=cumsum(w),1);
end
end
    
