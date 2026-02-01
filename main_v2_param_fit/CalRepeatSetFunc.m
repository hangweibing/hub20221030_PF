function [repeatSet] = CalRepeatSetFunc(spectra)
%CALREPEATSETFUNC 此处显示有关此函数的摘要
%   此处显示详细说明
if(mod(spectra,1)~=0)
    err('循环数不是整数')
end
cycleSet=spectra(1:2:end-1)-spectra(2:2:end);
tmp=1;
k=1;
for i=2:length(cycleSet)
    if (cycleSet(i)==cycleSet(i-1))
        tmp=tmp+1;
    else
        repeatSet(k)=tmp;
        tmp=1;
        k=k+1;
    end
end
end

