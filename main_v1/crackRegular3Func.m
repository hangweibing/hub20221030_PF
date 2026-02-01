function [yRegSet,zRegSet,ksiIniSet] = crackRegular3Func(yIniSet,zIniSet,nRegPoint,indexNeedIniNatCoord)
%CRACKREGULARFUNC 裂纹归一化程序，将非等间距非等结点的裂纹前缘归一化为等间距等结点裂纹
%   此处显示详细说明
% 下一步将列表改为元胞
% 这里需要确保y和z都是单调的，如果不单调，将会放大变形幅度
nDiscrPoint=10*nRegPoint;%建议样条曲线采样点比归一化的点多10倍以上，以避免在一些梯度大的地方出现两个采样点之间存在多个自然坐标点
% cs = spline(cfrNodeCoordList(:,2),cfrNodeCoordList(:,3));%样条函数
nCfrNode=size(yIniSet,1);
% zIniSet=cfrNodeCoordList(:,3)';
% yIniSet=cfrNodeCoordList(:,2)';
% cs = csape(yIniSet,zIniSet,'second');%样条函数
cs = csape(yIniSet,zIniSet);%样条函数
ySplineSet = linspace(yIniSet(1),yIniSet(end),nDiscrPoint);%插值点
zSplineSet=ppval(cs,ySplineSet);%插值
% plot(zIniSet,yIniSet,'r',zSplineSet,ySplineSet,'g'); 

%% 给出等间距地自然坐标，并给出原有点地自然坐标近似，这里的基础是ksi,y和z是一一对应的
nIntval=nRegPoint-1;
s_tot=0;
for i=2:nDiscrPoint
    ds=sqrt((ySplineSet(i)-ySplineSet(i-1)).^2+(zSplineSet(i)-zSplineSet(i-1)).^2);
    s_tot=s_tot+ds;
end

s_cur=0;
s_prv=0;
indexRegPoint=1;
indexIniPoint=2;
regArcSet=zeros(1,nRegPoint);
yRegSet=zeros(1,nRegPoint);
yRegSet(1)=yIniSet(1);
zRegSet=zeros(1,nRegPoint);
zRegSet(1)=zIniSet(1);
ksiIniSet=zeros(1,nCfrNode);
if(indexNeedIniNatCoord)
    for i=2:nDiscrPoint
        ds=sqrt((ySplineSet(i)-ySplineSet(i-1)).^2+(zSplineSet(i)-zSplineSet(i-1)).^2);
        s_prv=s_cur;
        s_cur=s_cur+ds;
        prvNatureCoord=s_prv/s_tot;
        curNatureCoord=s_cur/s_tot;
        reqNatureCoord=indexRegPoint/nIntval;
        % while语句，当出现存在两个样条曲线采样点之间存在两个或以上的自然坐标点，index+1，对新的自然坐标继续判断，并继续插值√
        while((prvNatureCoord<reqNatureCoord)&&(curNatureCoord>=reqNatureCoord))
            % 在两个点之间加线性插值，直接给出点的坐标
            yRegSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[ySplineSet(i-1),ySplineSet(i)],reqNatureCoord);
            zRegSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[zSplineSet(i-1),zSplineSet(i)],reqNatureCoord);
            regArcSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[s_prv,s_cur],reqNatureCoord)-sum(regArcSet(1:indexRegPoint));
            indexRegPoint=indexRegPoint+1;
            reqNatureCoord=indexRegPoint/nIntval;
        end
        % 根据y和z的坐标来确定自然参数
        while(indexIniPoint<=nCfrNode&&(ySplineSet(i-1)<yIniSet(indexIniPoint))&&(ySplineSet(i)>=yIniSet(indexIniPoint)))
            % 如果y和z均具有单值性，也就是说一个ksi对应一个y和一个z，那么使用y或者z寻找初始点对应的自然坐标都可以，恰巧这个问题的都满足
            %这里使用y
            ksiIniSet(indexIniPoint)=interp1([ySplineSet(i-1),ySplineSet(i)],[prvNatureCoord,curNatureCoord],yIniSet(indexIniPoint));
            indexIniPoint=indexIniPoint+1;
        end
    end
else
    for i=2:nDiscrPoint
        ds=sqrt((ySplineSet(i)-ySplineSet(i-1)).^2+(zSplineSet(i)-zSplineSet(i-1)).^2);
        s_prv=s_cur;
        s_cur=s_cur+ds;
        prvNatureCoord=s_prv/s_tot;
        curNatureCoord=s_cur/s_tot;
        reqNatureCoord=indexRegPoint/nIntval;
        % while语句，当出现存在两个样条曲线采样点之间存在两个或以上的自然坐标点，index+1，对新的自然坐标继续判断，并继续插值√
        while((prvNatureCoord<reqNatureCoord)&&(curNatureCoord>=reqNatureCoord))
            % 在两个点之间加线性插值，直接给出点的坐标
            yRegSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[ySplineSet(i-1),ySplineSet(i)],reqNatureCoord);
            zRegSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[zSplineSet(i-1),zSplineSet(i)],reqNatureCoord);
            regArcSet(indexRegPoint+1)=interp1([prvNatureCoord,curNatureCoord],[s_prv,s_cur],reqNatureCoord)-sum(regArcSet(1:indexRegPoint));
            indexRegPoint=indexRegPoint+1;
            reqNatureCoord=indexRegPoint/nIntval;
        end
    end
end
if(zRegSet(end)<=35&&yRegSet(1)>=7)
    ksiRegSet=linspace(0,1,nRegPoint);
    yz=[yRegSet;zRegSet];
    cs5=spline(ksiRegSet,[[0;1],yz,[1;0]]);
    yzNew=ppval(cs5,linspace(ksiRegSet(1),ksiRegSet(end),nRegPoint));
    % figure(5)
    % plot(XY(1,:),XY(2,:),'bo',yy(1,:),yy(2,:))
    yRegSet=yzNew(1,:);
    zRegSet=yzNew(2,:);
end
end

