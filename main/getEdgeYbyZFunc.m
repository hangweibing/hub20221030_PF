function [y] = getEdgeYbyZFunc(z,location)
%GETEDGEYBYZFUNC 此处显示有关此函数的摘要
%   此处显示详细说明
if(strcmp(location,'up'))
    if(z<=35)
        y=13;
    end
    if(z>35&&z<=37.82842712)
        y=14-sqrt(9-(z-37.82842712).^2);
    end
    if(z>37.82842712)
        y=11;
    end
end
if(strcmp(location,'down'))
    if(z<=35)
        y=7; 
    end
    if(z>35&&z<=37.82842712)
        y=sqrt(9-(z-37.82842712).^2)+6;
    end
    if(z>37.82842712)
        y=9;
    end
    if(~exist('y','var'))
        y=7;
    end
        
end
end

