xparticleForCalII=xparticle;
% ac=32.5; %用重采样后的粒子计算检查间隔，保留材料参数不变
a_ini_set_forIICal=(wblrnd(0.0057,1.1739,1,N)+0.0446369).*25.4;
for i=1:N
    a=a_ini_set_forIICal(i);
    c=a;
    for j=1:n_nodes
        yIniRegSet(j)=centers(1)+c*sin(thetas(j));
        zIniRegSet(j)=centers(2)+a*cos(thetas(j));
    end
    xparticleForCalII(i,1:2*n_nodes)=[yIniRegSet,zIniRegSet];
end
SPLITTE_temp=zeros(1,N);
xparticlek_1=xparticleForCalII;
aver_delta_sigma=aver_delta_sigma_set(1);
parfor i=1:N
    curUinput={};
    curAverInput={};
    SPLITTED=SPLITTE_temp(i);
    xparticlei=xparticlek_1(i,:);
    a_up=xparticlei(42)-30;
    a_down=xparticlei(22)-30;
    m_index=0;
    stage=[1,2,3,4,5,6,7,0,0,0,0,8];% 每个阶段在元胞中的位置
    while m_index==0
        m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
        if m_index==0
            tmp_sel=randi([1,N],1,1);
            xparticlei=xparticlek_1(tmp_sel,:);
            a_up=xparticlei(42)-30;
            a_down=xparticlei(22)-30;
        end
    end
    if SPLITTED && (m_index==3 || m_index==5)  % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
        if m_index==3
            curUinput=Uinput_splitted_1;
            curAverInput=averInput_splitted_1;
        elseif m_index==5
            curUinput=Uinput_splitted_2;
            curAverInput=averInput_splitted_2;
        end
        m_name=sprintf('nn_stage%ds',stage(m_index));
    else % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
        curUinput=Uinput_integrated{stage(m_index)};
        curAverInput=averInput_integrated{stage(m_index)};
        m_name=sprintf('nn_stage%d',stage(m_index));
    end
    yRegSet=xparticlei(1:21);
    zRegSet=xparticlei(22:42);
    logCstar=xparticlei(43);
    gamma=xparticlei(44);
    [yRegSet,zRegSet,SPLITTED,logCstar,gamma,deltaK] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
    deltaKSet(i)=deltaK;
end
[y1,x1]=ksdensity(deltaKSet);
[y1,x1]=ksdensity(deltaKSet,linspace(min(x1),max(x1)+3,N));
label=find(x1>0);
xs=x1(label);
ys=y1(label);
yr=normcdf(xs,33.4,3.34);
yrs=yr.*ys;
pofCum=cumtrapz(xs,yrs);
pof=pofCum(end);
PoFSet(end)=pof;

% [y1,x1]=ksdensity(xparticleForCalII(:,42));
% [y1,x1]=ksdensity(xparticleForCalII(:,42),linspace(min(x1)-0.5,max(x1)+0.5,N));  % N为粒子数目
% csap=csape(x1,y1);
% if ac<max(x1)
%     pdf_ac=ppval(csap,ac);
%     n1=sum(x1<ac);
%     x2=[ac x1(n1+1:end)];
%     y2=[pdf_ac y1(n1+1:end)];
%     pofCum=cumtrapz(x2,y2);
%     pof=pofCum(end);
% else
%     pof=0;
% end
% PoFSet(end)=pof;

%%
SPLITTE_temp=zeros(1,N);
j=2;
k1=j-1;
% deltaKSet=zeros(1,N);
while PoFSet(end)<10^(-7)
    xparticlek_1=xparticleForCalII;
    aver_delta_sigma=aver_delta_sigma_set(j-1);
    parfor i=1:N
        curUinput={};
        curAverInput={};
        SPLITTED=SPLITTE_temp(i);
        xparticlei=xparticlek_1(i,:);
        a_up=xparticlei(42)-30;
        a_down=xparticlei(22)-30;
        m_index=0;
        stage=[1,2,3,4,5,6,7,0,0,0,0,8];% 每个阶段在元胞中的位置
        while m_index==0
            m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
            if m_index==0
                tmp_sel=randi([1,N],1,1);
                xparticlei=xparticlek_1(tmp_sel,:);
                a_up=xparticlei(42)-30;
                a_down=xparticlei(22)-30;
            end
        end
        if SPLITTED && (m_index==3 || m_index==5)  % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
            if m_index==3
                curUinput=Uinput_splitted_1;
                curAverInput=averInput_splitted_1;
            elseif m_index==5
                curUinput=Uinput_splitted_2;
                curAverInput=averInput_splitted_2;
            end
            m_name=sprintf('nn_stage%ds',stage(m_index));
        else % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
            curUinput=Uinput_integrated{stage(m_index)};
            curAverInput=averInput_integrated{stage(m_index)};
            m_name=sprintf('nn_stage%d',stage(m_index));
        end
        yRegSet=xparticlei(1:21);
        zRegSet=xparticlei(22:42);
        logCstar=xparticlei(43);
        gamma=xparticlei(44);
        [yRegSet,zRegSet,SPLITTED,logCstar,gamma,deltaK] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
        xparticleForCalII(i,:)=real([yRegSet,zRegSet,logCstar,gamma]);
        SPLITTE_temp(i)=SPLITTED;
        deltaKSet(i)=deltaK;
    end
    [y1,x1]=ksdensity(deltaKSet);
    [y1,x1]=ksdensity(deltaKSet,linspace(min(x1),max(x1)+3,N));
    label=find(x1>0);
    xs=x1(label);
    ys=y1(label);
    yr=normcdf(xs,33.4,3.34);
    yrs=yr.*ys;
    pofCum=cumtrapz(xs,yrs);
    pof=pofCum(end);
    %     [y1,x1]=ksdensity(xparticleForCalII(:,42),linspace(min(x1)-0.5,max(x1)+0.5,N));  % N为粒子数目
    %     csap=csape(x1,y1);
    %     if ac<max(x1)
    %         pdf_ac=ppval(csap,ac);
    %         n1=sum(x1<ac);
    %         x2=[ac x1(n1+1:end)];
    %         y2=[pdf_ac y1(n1+1:end)];
    %         pofCum=cumtrapz(x2,y2);
    %         pof=pofCum(end);
    %     else
    %         pof=0;
    %     end
    PoFSet=[PoFSet pof];
    j=j+1;
    disp(['检查间隔计算：第' num2str((j-1)*step/1950.70866) '个小时; 失效概率为' num2str(PoFSet(end))]);
end
k2=j-1;
InspectionInterval=(k2-k1)*step/1950.70866;