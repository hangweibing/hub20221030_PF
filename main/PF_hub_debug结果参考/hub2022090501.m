%% preliminary treatment
clear;close all;clc;tic;format long;
% addpath surModelPackage
% 正则粒子滤波初始化
n=1;                                    %状态向量中每一个随机变量的维数
N=1;                                 %粒子数目
v_sphere=2;                             %一维空间球体体积
A=(8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h=A*N^(-1/(n+4));                       %参见文献《Robust regularized particle filter for terraino.k;  vvhjk ,. navigation》

D=zeros(44,1000);                       %经验方差的开根值
e=zeros(N,44,1000);                     %从Epanechikov核函数中的采样结果
cyclesperhour=1950.70866;
Xpf=zeros(44,1000);                     %每一个时间步的估计值
xparticle=zeros(N,44,1000);
xparticle1=zeros(N,44,1000);
xparticle_cov=zeros(44,44,1000);
upcrackparticles=zeros(N,1000);
logCstarparticles=zeros(N,1000);
gammaparticles=zeros(N,1000);
weight=zeros(N,1000);
% t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
% z      =[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;
t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
z      =[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;
zPred=zeros(N,1);
R=0.5;
n_nodes=21;
centers=[13,30];
downCenter=[6,37.82842712]; % 下圆弧圆心
downRadius=3; % 下圆弧半径
nRegPoint=n_nodes; % 等间距点数目
thetas=linspace(3/2*pi,2*pi,n_nodes);
yIniRegSet=zeros(1,n_nodes);
zIniRegSet=zeros(1,n_nodes);
for i=1:N   %粒子集初始化
%     logCstar=unifrnd(-10.9,-10.8,1,1);
    logCstar=-10.8;
%     gamma=normrnd(3,0.01,1,1);
    gamma=3.0;
%     a=normrnd(2,0.01,1,1);
    a=0.1;
    c=a;
    for j=1:n_nodes
        yIniRegSet(j)=centers(1)+c*sin(thetas(j));
        zIniRegSet(j)=centers(2)+a*cos(thetas(j));
    end
    xparticle(i,:,1)=[yIniRegSet,zIniRegSet,logCstar,gamma];
end
a_up=a;
upcrackparticles(:,1)=xparticle(:,42,1);
logCstarparticles(:,1)=xparticle(:,43,1);
gammaparticles(:,1)=xparticle(:,44,1);
weight(:,1)=1/N*ones(N,1);
Xpf(:,1)=(mean(xparticle(:,:,1)))';
xparticle_cov(:,:,1)=cov(xparticle(:,:,1));

load('AsteixSpectraData.mat');
load('pod_models.mat'); %加载代理模型数据
spectra=[spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra];
spectra=spectra(2:end);
step=1000;
aver_delta_sigma_set=zeros(1,ceil((length(spectra))/2/step));
k=1;
for i=1:floor((length(spectra))/2/step)
    delta_sigmas=zeros(1,step);
    for j=1:step
        Smax=spectra(2*(k+j-1));
        Smin=spectra(2*(k+j-1)-1);
        delta_sigmas(j)=Smax-Smin;
    end
    aver_delta_sigma_set(i)=mean(delta_sigmas);
    k=k+step;
end
%


delta_sigmas=0;
k=floor((length(spectra))/2/step)*step+1;
i=1;
while k<=(length(spectra)/2)
    Smax=spectra(2*k);
    Smin=spectra(2*k-1);
    delta_sigmas(i)=Smax-Smin;
    k=k+1;
    i=i+1;
end
aver_delta_sigma_set(end)=mean(delta_sigmas);
stage=[1,2,3,4,5,6,7,0,0,0,0,8];
testErrSet=[0.019803420755871 0.058377623733943 0.013343080193008 0.009221345729625 0.037283991994545 0.011408674471790 0.082711915490345 0.032375406062069 0.041104885753232 0.057544490601307];
%
clc
plot_geometry_20
hold on
SPLITTE_temp=zeros(1,N);
m=2;j=1;
% while (m-1)*step/1950.70866<=t_check(end)
while a_up<19.8
    xparticlem_1=xparticle(:,:,m-1);
    aver_delta_sigma=aver_delta_sigma_set(m-1);
    for i=1:N
        curUinput={};
        curAverInput={};
        SPLITTED=SPLITTE_temp(i);
        xparticlei=xparticlem_1(i,:);
        a_up=xparticlei(42)-30;
        a_down=xparticlei(22)-30;
        m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
        if ~SPLITTED % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
            curUinput=Uinput_integrated{stage(m_index)};
            curAverInput=averInput_integrated{stage(m_index)};
            m_name=sprintf('nn_stage%d',stage(m_index));
        else % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
            if m_index==3
                curUinput=Uinput_splitted{1};
                curAverInput=averInput_splitted{1};
            elseif m_index==5
                curUinput=Uinput_splitted{2};
                curAverInput=averInput_splitted{2};
            end
            m_name=sprintf('nn_stage%ds',stage(m_index));
        end
        yRegSet=xparticlei(1:21);
        zRegSet=xparticlei(22:42);
        logCstar=xparticlei(43);
        gamma=xparticlei(44);
        [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aFunc(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
        if mod(m,2)==0
            plot(zRegSet,yRegSet);
            hold on
        end
        a_up=zRegSet(end)-30;
        fprintf('a_up: %f a_down: ',zRegSet(end)-30,zRegSet(1)-30)
        xparticle(i,:,m)=[yRegSet,zRegSet,logCstar,gamma];
        SPLITTE_temp(i)=SPLITTED;
    end
    weight(:,m)=1/N*ones(N,1);
    Xpf(:,m)=(mean(xparticle(:,:,m)))';
    xparticle_cov(:,:,m)=cov(xparticle(:,:,m));
    upcrackparticles(:,m)=xparticle(:,42,m);
    logCstarparticles(:,m)=xparticle(:,43,m);
    gammaparticles(:,m)=xparticle(:,44,m);
%     if j<=4
%         R=0.5
%     else
%         R=0.3
%     end
%     if ((m-1)*step/1950.70866<t_check(j))&&(m*step/1950.70866>=t_check(j))
%         parfor i=1:N   %向后推理
%             zPred(i)=upcrackparticles(i,m);
%             z1(i)=z(j)-zPred(i);
%             weight(i,m)=inv(sqrt(2*pi*det(R)))*exp(-0.5*(z1(i))*inv(R)*(z1(i))')+1e-99;
%         end
%         weight(:,m)=weight(:,m)./sum(weight(:,m));
%         Xpf(:,m)=0;
%         for i=1:N
%             Xpf(:,m)=Xpf(:,m)+(weight(i,m)*xparticle(i,:,m))';
%         end
%         xparticle_cov(:,:,m)=0;
%         for i=1:N
%             xparticle_cov(:,:,m)=xparticle_cov(:,:,m)+weight(i,m)*(xparticle(i,:,m)'-Xpf(:,m))*(xparticle(i,:,m)'-Xpf(:,m))';
%         end
%         for i=1:44
%             D(i,m)=sqrt(xparticle_cov(i,i,m));
%             e(:,i,m)=kernelsampling(N)';
%         end
%         outindex=randomr(weight(:,m));
%         xparticle1(:,:,m)=xparticle(outindex,:,m);                                 %重采样
%         for i=1:44
%             xparticle(:,i,m)=xparticle1(:,i,m)+h*D(i,m)*e(:,i,m);              %正则化
%             xparticle(:,i,m)=rearrange(xparticle(:,i,m),xparticle1(:,i,m))';     %参见文献《Dynamic Bayesian Network for Aircraft Wing Health Monitoring Digital Twin》
%         end
%         for i=1:N
%             [yNewSet,zNewSet,SPLITTE_temp(i)]=addConstraintNewSatgeFunc(xparticle(i,1:21,m),xparticle(i,22:42,m));
%             [xparticle(i,1:21,m),xparticle(i,22:42,m),~] = crackRegular5Func(yNewSet,zNewSet,nRegPoint,'false');
%         end
%         j=j+1;
%     end
    m=m+1;
    disp(['第' num2str((m-1)*step/1950.70866) '个小时，经历了' num2str(j-1) '次检查']);
end
toc
%%
close all
figure
clear x y_0 y_1 y_2 y_3 y_33 y21 PoF3
for i=1:m-1
    y_1(i)=prctile(upcrackparticles(:,i),99.95)-30;
    y_2(i)=prctile(upcrackparticles(:,i),0.05)-30;
    y_3(i)=Xpf(42,i)-30;
    x(i)=(i-1)*step/1950.70866;
end
t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
z=[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;

% plot(t_check,z-30,'^','linewidth',5);hold on
% plot(x,y_3,'b','linewidth',5);hold on;
% plot(x,y_1,'r--','linewidth',5);hold on
plot(x,y_2,'r--','linewidth',5);hold off;


xlabel('Flight hours/h','FontSize',30);
ylabel('Surface crack length/mm','FontSize',30);
set(get(gca,'xlabel'),'fontname','Times New Roman');
set(get(gca,'ylabel'),'fontname','Times New Roman');
set(gca,'fontname','Times New Roman');
set(gca,'FontSize',40);
% legend('   Experimental value','   Prediction mean', '   99.9% bounds','FontSize',40);
legend('boxoff')
legend('Location','best');
grid on;
set(gca,'gridlinestyle',':','gridcolor','k');