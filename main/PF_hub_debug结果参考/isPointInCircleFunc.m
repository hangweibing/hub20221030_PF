function judge=isPointInCircleFunc(point,center,r)
% 判断点是否在圆内
distance=sqrt((point(1)-center(1)).^2+(point(2)-center(2)).^2);
if(distance<=r)
    judge=true;
else
    judge=false;
end
end
