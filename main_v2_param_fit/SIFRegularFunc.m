function [K1RegSet,outputArg2] = SIFRegularFunc(inputArg1,inputArg2)
%SIFREGULARFUNC 此处显示有关此函数的摘要
%   此处显示详细说明
csk1Ini = csape(ksiIniSet,K1Set,'second');%样条函数
csk2Ini = csape(ksiIniSet,K2Set,'second');%样条函数
K1RegSet=ppval(csk1Ini,ksiRegSet);%插值
K2RegSet=ppval(csk2Ini,ksiRegSet);%插值
end

