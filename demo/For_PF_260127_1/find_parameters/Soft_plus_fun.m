function [rou] = Soft_plus_fun(alpha,x)
%SOFT_PLUS_FUN 此处显示有关此函数的摘要

    rou=1/alpha*log( 1+exp(alpha*x) ); 


end

