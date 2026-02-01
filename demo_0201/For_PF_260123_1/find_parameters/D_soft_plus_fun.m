function [D_rou] = D_soft_plus_fun(alpha,x)

    D_rou=1./( 1+exp(-1*alpha*x) );

end

