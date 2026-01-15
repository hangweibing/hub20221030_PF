function [yNewSet,zNewSet,splitted]=resplit_front(yNewSet,zNewSet,downCenter,r)
%当裂纹前缘与下圆角相交时，给出裂纹前缘与几何边界的合并曲线

splitted=false;
n_reg_point=length(yNewSet);
judge=zeros(1,n_reg_point);
for i=1:length(yNewSet)
    judge(i)=isPointInCircleFunc([yNewSet(i),zNewSet(i)],downCenter,r);
end
intersected_coords=find(judge==true);
if (~isempty(intersected_coords))
    splitted=true;
    left_coord=intersected_coords(1);
    right_coord=intersected_coords(end);
    if left_coord~=1
        [leftInsecY,leftInsecZ]=getPointOnCurveFunc([yNewSet(left_coord-1),zNewSet(left_coord-1)]...
            ,[yNewSet(left_coord),zNewSet(left_coord)],downCenter,r,'left');
        [rightInsecY,rightInsecZ]=getPointOnCurveFunc([yNewSet(right_coord),zNewSet(right_coord)]...
        ,[yNewSet(right_coord+1),zNewSet(right_coord+1)],downCenter,r,'right');
        y_add=linspace(leftInsecY,rightInsecY,10);
        z_add=linspace(leftInsecZ,rightInsecZ,10);

        yNewSet=[yNewSet(1:left_coord-1),y_add,yNewSet(right_coord+1:end)];
        zNewSet=[zNewSet(1:left_coord-1),z_add,zNewSet(right_coord+1:end)];
    else
        % 当下断点已经不在下边界时，只取上半段，这是需要对结点进行修正,但是下端点不一定在下圆弧上，有可能已经到右下边界了
        if zNewSet(right_coord)<=37.82842712
            [rightInsecY,rightInsecZ]=getPointOnCurveFunc([yNewSet(right_coord),zNewSet(right_coord)]...
            ,[yNewSet(right_coord+1),zNewSet(right_coord+1)],downCenter,r,'right');
        else
            rightInsecZ=interp1([yNewSet(right_coord),yNewSet(right_coord+1)],[zNewSet(right_coord),zNewSet(right_coord+1)],9,'linear','extrap');%也有可能出现外插
            rightInsecY=9;
        end
        splitted=false;
        yNewSet=[rightInsecY,yNewSet(right_coord+1:end)];
        zNewSet=[rightInsecZ,zNewSet(right_coord+1:end)];
    end    
end
end
