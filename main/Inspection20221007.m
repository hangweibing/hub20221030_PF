%% preliminary treatment
clear all
close all
clc
tic
load('AsteixSpectraData.mat');
load('pod_models.mat'); %加载代理模型数据
n_nodes=21;
centers=[13,30];
downCenter=[6,37.82842712]; % 下圆弧圆心
downRadius=3; % 下圆弧半径
nRegPoint=n_nodes; % 等间距点数目
thetas=linspace(3/2*pi,2*pi,n_nodes);
yIniRegSet=zeros(1,n_nodes);
zIniRegSet=zeros(1,n_nodes);
a_ini_set=lognrnd(-2.65,0.297,1,1000);
a_ini_test=a_ini_set(randi([1,1000],1,4));
mu=[-10.9 3];
SIGMA=[0.15^2 -0.05^2;-0.05^2 0.05^2];
paraSet=mvnrnd(mu,SIGMA,1000);

%
% 试件1初始化
test_coupon_1=[a_ini_test(1) paraSet(randi([1,1000],1,1),:)];
xparticle_coupon1=zeros(1,44);
a=test_coupon_1(1);
c=a;
for j=1:n_nodes
    yIniRegSet(j)=centers(1)+c*sin(thetas(j));
    zIniRegSet(j)=centers(2)+a*cos(thetas(j));
end
xparticle_coupon1(1,:)=[yIniRegSet,zIniRegSet,test_coupon_1(2),test_coupon_1(3)];

% 试件2初始化
test_coupon_2=[a_ini_test(2) paraSet(randi([1,1000],1,1),:)];
xparticle_coupon2=zeros(1,44);
a=test_coupon_2(1);
c=a;
for j=1:n_nodes
    yIniRegSet(j)=centers(1)+c*sin(thetas(j));
    zIniRegSet(j)=centers(2)+a*cos(thetas(j));
end
xparticle_coupon2(1,:)=[yIniRegSet,zIniRegSet,test_coupon_2(2),test_coupon_2(3)];

% 试件3初始化
test_coupon_3=[a_ini_test(3) paraSet(randi([1,1000],1,1),:)];
xparticle_coupon3=zeros(1,44);
a=test_coupon_3(1);
c=a;
for j=1:n_nodes
    yIniRegSet(j)=centers(1)+c*sin(thetas(j));
    zIniRegSet(j)=centers(2)+a*cos(thetas(j));
end
xparticle_coupon3(1,:)=[yIniRegSet,zIniRegSet,test_coupon_3(2),test_coupon_3(3)];

% 试件4初始化
test_coupon_4=[a_ini_test(4) paraSet(randi([1,1000],1,1),:)];
xparticle_coupon4=zeros(1,44);
a=test_coupon_4(1);
c=a;
for j=1:n_nodes
    yIniRegSet(j)=centers(1)+c*sin(thetas(j));
    zIniRegSet(j)=centers(2)+a*cos(thetas(j));
end
xparticle_coupon4(1,:)=[yIniRegSet,zIniRegSet,test_coupon_4(2),test_coupon_4(3)];

spectra=[spectra spectra spectra spectra spectra spectra spectra spectra spectra spectra];
spectra=spectra(2:end);
step=2000;
aver_delta_sigma_set=zeros(1,ceil((length(spectra))/2/step));
k=1;  % 载荷谱里的载荷循环数
for i=1:floor((length(spectra))/2/step)
    delta_sigmas=zeros(1,step);
    for j =1:step
        Smax=spectra(2*(k+j-1));
        Smin=spectra(2*(k+j-1)-1);
        delta_sigmas(j)=Smax-Smin;
    end
    aver_delta_sigma_set(i)=mean(delta_sigmas);
    k=k+step;
end
ac=33.44;
Uinput_splitted_1=Uinput_splitted{1};
averInput_splitted_1=averInput_splitted{1};
Uinput_splitted_2=Uinput_splitted{2};
averInput_splitted_2=averInput_splitted{2};
testErrSet=[0.019803420755871 0.058377623733943 0.013343080193008 0.009221345729625 0.037283991994545 0.011408674471790 0.082711915490345 0.032375406062069 0.041104885753232 0.057544490601307];
%
% 计算第一个试件的真实裂纹扩展历程
SPLITTE_temp=0;
i1=1;
while xparticle_coupon1(i1,42)<ac
    disp(['试件1，裂纹长度=' num2str(xparticle_coupon1(i1,42))]);
    aver_delta_sigma=aver_delta_sigma_set(i1);
    curUinput={};
    curAverInput={};
    SPLITTED=SPLITTE_temp;
    a_up=xparticle_coupon1(i1,42)-30;
    a_down=xparticle_coupon1(i1,22)-30;
    m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
    stage=[1,2,3,4,5,6,7,0,0,0,0,8]; % 每个阶段在元胞中的位置
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
    yRegSet=xparticle_coupon1(i1,1:21);
    zRegSet=xparticle_coupon1(i1,22:42);
    logCstar=xparticle_coupon1(i1,43);
    gamma=xparticle_coupon1(i1,44);
    [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
    i1=i1+1;
    xparticle_coupon1(i1,:)=[yRegSet,zRegSet,logCstar,gamma];
    SPLITTE_temp=SPLITTED;
end
%
% 计算第二个试件的真实裂纹扩展历程
SPLITTE_temp=0;
i2=1;
while xparticle_coupon2(i2,42)<ac
    disp(['试件2，裂纹长度=' num2str(xparticle_coupon2(i2,42))]);
    aver_delta_sigma=aver_delta_sigma_set(i2);
    curUinput={};
    curAverInput={};
    SPLITTED=SPLITTE_temp;
    a_up=xparticle_coupon2(i2,42)-30;
    a_down=xparticle_coupon2(i2,22)-30;
    m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
    stage=[1,2,3,4,5,6,7,0,0,0,0,8]; % 每个阶段在元胞中的位置
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
    yRegSet=xparticle_coupon2(i2,1:21);
    zRegSet=xparticle_coupon2(i2,22:42);
    logCstar=xparticle_coupon2(i2,43);
    gamma=xparticle_coupon2(i2,44);
    [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
    i2=i2+1;
    xparticle_coupon2(i2,:)=[yRegSet,zRegSet,logCstar,gamma];
    SPLITTE_temp=SPLITTED;
end

% 计算第三个试件的真实裂纹扩展历程
SPLITTE_temp=0;
i3=1;
while xparticle_coupon3(i3,42)<ac
    disp(['试件3，裂纹长度=' num2str(xparticle_coupon3(i3,42))]);
    aver_delta_sigma=aver_delta_sigma_set(i3);
    curUinput={};
    curAverInput={};
    SPLITTED=SPLITTE_temp;
    a_up=xparticle_coupon3(i3,42)-30;
    a_down=xparticle_coupon3(i3,22)-30;
    m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
    stage=[1,2,3,4,5,6,7,0,0,0,0,8]; % 每个阶段在元胞中的位置
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
    yRegSet=xparticle_coupon3(i3,1:21);
    zRegSet=xparticle_coupon3(i3,22:42);
    logCstar=xparticle_coupon3(i3,43);
    gamma=xparticle_coupon3(i3,44);
    [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
    i3=i3+1;
    xparticle_coupon3(i3,:)=[yRegSet,zRegSet,logCstar,gamma];
    SPLITTE_temp=SPLITTED;
end

% 计算第四个试件的真实裂纹扩展历程
SPLITTE_temp=0;
i4=1;
while xparticle_coupon4(i4,42)<ac
    disp(['试件4，裂纹长度=' num2str(xparticle_coupon4(i4,42))]);
    aver_delta_sigma=aver_delta_sigma_set(i4);
    curUinput={};
    curAverInput={};
    SPLITTED=SPLITTE_temp;
    a_up=xparticle_coupon4(i4,42)-30;
    a_down=xparticle_coupon4(i4,22)-30;
    m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
    stage=[1,2,3,4,5,6,7,0,0,0,0,8]; % 每个阶段在元胞中的位置
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
    yRegSet=xparticle_coupon4(i4,1:21);
    zRegSet=xparticle_coupon4(i4,22:42);
    logCstar=xparticle_coupon4(i4,43);
    gamma=xparticle_coupon4(i4,44);
    [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aNew(yRegSet,zRegSet,aver_delta_sigma,m_name,curUinput,curAverInput,logCstar,gamma,step,testErrSet);
    i4=i4+1;
    xparticle_coupon4(i4,:)=[yRegSet,zRegSet,logCstar,gamma];
    SPLITTE_temp=SPLITTED;
end
%%
clear;load couponBase20220921.mat
figure
x=zeros(1,i1);
y=zeros(1,i1);
for i=1:i1
    x(i)=(i-1)*step/1950.70866;
    y(i)=xparticle_coupon1(i,42)-30;
end
n1=length(y(y<1.49));n2=length(y(y<2.9));
cs=csape([x(1:n1) x(n2:end)],[y(1:n1) y(n2:end)]);
y1=ppval(cs,x);coupon1TimeStep=x;coupon1ActualCrack=y1;
% plot(x,y,'linewidth',2);hold on
plot(x,y1,'linewidth',2);hold on
% legend('before','after')
%
x=zeros(1,i2);
y=zeros(1,i2);
for i=1:i2
    x(i)=(i-1)*step/1950.70866;
    y(i)=xparticle_coupon2(i,42)-30;
end
n1=length(y(y<1.49));n2=length(y(y<2.9));
cs=csape([x(1:n1) x(n2:end)],[y(1:n1) y(n2:end)]);
y1=ppval(cs,x);coupon2TimeStep=x;coupon2ActualCrack=y1;
plot(x,y1,'linewidth',2);hold on


%
x=zeros(1,i3);
y=zeros(1,i3);
for i=1:i3
    x(i)=(i-1)*step/1950.70866;
    y(i)=xparticle_coupon3(i,42)-30;
end
n1=length(y(y<1.49));n2=length(y(y<2.9));
cs=csape([x(1:n1) x(n2:end)],[y(1:n1) y(n2:end)]);
y1=ppval(cs,x);coupon3TimeStep=x;coupon3ActualCrack=y1;
plot(x,y1,'linewidth',2);hold on


%
x=zeros(1,i4);
y=zeros(1,i4);
for i=1:i4
    x(i)=(i-1)*step/1950.70866;
    y(i)=xparticle_coupon4(i,42)-30;
end
n1=length(y(y<1.49));n2=length(y(y<2.9));
cs=csape([x(1:n1) x(n2:end)],[y(1:n1) y(n2:end)]);
y1=ppval(cs,x);coupon4TimeStep=x;coupon4ActualCrack=y1;
plot(x,y1,'linewidth',2);hold on

x=zeros(1,i6);
y=zeros(1,i6);
for i=1:i6
    x(i)=(i-1)*step/1950.70866;
    y(i)=coupon6ActualCrack(i);
end
% n1=length(y(y<1.49));n2=length(y(y<2.9));
% cs=csape([x(1:n1) x(n2:end)],[y(1:n1) y(n2:end)]);
% y1=ppval(cs,x);coupon1TimeStep=x;coupon1ActualCrack=y1;
% plot(x,y,'linewidth',2);hold on
plot(x,y,'linewidth',2);hold on


hold off
xlabel('飞行时间/h','FontSize',20);
ylabel('裂纹长度/mm','FontSize',20);
set(gca,'FontSize',20);
legend('试件1','试件2','试件3','试件4','试件6','FontSize',20);
legend('Location','North');

%% 正则粒子滤波初始化
clear;clc;
n=1;                                    %状态向量中每一个随机变量的维数
N=200000;                                %粒子数目
v_sphere=2;                             %一维空间球体体积
A=(8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h=A*N^(-1/(n+4));                       %参见文献《Robust regularized particle filter for terraino.k;  vvhjk ,. navigation》
D=zeros(44,3000);                       %经验方差的开根值
e=zeros(N,44,1);                        %从Epanechikov核函数中的采样结果
cyclesperhour=1950.70866;
Xpf=zeros(44,3000);                     %每一个时间步的估计值
xparticle=zeros(N,44,1);
xparticle1=zeros(N,44,1);
xparticle_cov=zeros(44,44,1);
upcrackparticles=zeros(N,3000);
logCstarparticles=zeros(N,3000);
gammaparticles=zeros(N,3000);
weight=zeros(N,3000);
zPred=zeros(N,1);
R=0.01;
n_nodes=21;
yIniRegSet=zeros(1,n_nodes);
zIniRegSet=zeros(1,n_nodes);
%
load('AsteixSpectraData.mat');
load('pod_models.mat'); %加载代理模型数据
n_nodes=21;
centers=[13,30];
downCenter=[6,37.82842712]; % 下圆弧圆心
downRadius=3; % 下圆弧半径
nRegPoint=n_nodes; % 等间距点数目
thetas=linspace(3/2*pi,2*pi,n_nodes);

a_ini_Set=lognrnd(-2.65,0.297,1,N);
% a_ini_Set=normrnd(2,0.1,1,N);
centers=[13,30];
mu=[-10.9 3];
SIGMA=[0.15^2 -0.05^2;-0.05^2 0.05^2];
for i=1:N   %粒子集初始化
    paraSet=mvnrnd(mu,SIGMA,1);
    logCstar=paraSet(1);
    gamma=paraSet(2);
    a=a_ini_Set(i);
    c=a;
    for j=1:n_nodes
        yIniRegSet(j)=centers(1)+c*sin(thetas(j));
        zIniRegSet(j)=centers(2)+a*cos(thetas(j));
    end
    xparticle(i,:,1)=[yIniRegSet,zIniRegSet,logCstar,gamma];
end
upcrackparticles(:,1)=xparticle(:,42,1);
logCstarparticles(:,1)=xparticle(:,43,1);
gammaparticles(:,1)=xparticle(:,44,1);
weight(:,1)=1/N*ones(N,1);
Xpf(:,1)=(mean(xparticle(:,:,1)))';
xparticle_cov(:,:,1)=cov(xparticle(:,:,1));
SPLITTE_temp=zeros(1,N);

% load inspectionBasedata20220212 csap
% save inspectionBasedata20220223
load couponBase20220921.mat
SPLITTE_temp=zeros(1,N);
% %% firstInspection
%
tic

InspectionIntervalSet=[];PoFSet=0;
PoFBasedInspectionIntervalCal
InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
InspectionTime=sum(InspectionIntervalSet);
k=2;
SPLITTE_temp=zeros(1,N);
probabilisticPropagation
podX=[0.00255754	0.0383632	0.0869565	0.132992	0.171355	0.194373	0.214834	0.237852	0.258312	0.276215	0.29156	0.304348	0.319693	0.329923	0.340153	0.352941	0.365729	0.375959	0.383632	0.391304	0.396419	0.401535	0.40665	0.409207	0.41688	0.424552	0.429668	0.434783	0.442455	0.44757	0.455243	0.460358	0.465473	0.473146	0.480818	0.488491	0.498721	0.506394	0.516624	0.524297	0.531969	0.542199	0.55243	0.56266	0.578005	0.59335	0.606138	0.621483	0.636829	0.657289	0.677749	0.69821	0.721228	0.744246	0.764706	0.790281	0.83376	0.890026	0.933504	0.97954	1.04859	1.11509	1.19182	1.25831	1.30691	1.43734];
podY=[0	0	0	0	0	0.00451128	0.00902256	0.0180451	0.0315789	0.0451128	0.0616541	0.081203	0.105263	0.12782	0.15188	0.183459	0.218045	0.246617	0.273684	0.294737	0.314286	0.329323	0.347368	0.362406	0.381955	0.409023	0.431579	0.451128	0.47218	0.490226	0.515789	0.535338	0.550376	0.57594	0.601504	0.62406	0.652632	0.675188	0.697744	0.718797	0.736842	0.760902	0.780451	0.8	0.825564	0.84812	0.863158	0.879699	0.896241	0.911278	0.926316	0.938346	0.948872	0.957895	0.96391	0.971429	0.978947	0.986466	0.989474	0.992481	0.995489	0.996992	0.998496	0.998496	1	1];
csap=csape(podX,podY);
save firstInspectionTest20221007N_100000
toc
%% 试件1
% Intelligent inspection
% 试件1

clear;clc;load('firstInspectionTest20221001.mat');
%%
tic
inspectionValues=coupon1ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3;
trueValues=coupon1ActualCrack(k-1);
detectionNumber=10;

while inspectionValues(end)<1.27
    PFWithoutInspectionValues;
    PoFBasedInspectionIntervalCal;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon1ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon1ActualCrack(k-1)];
end
%

while inspectionValues(end)<2.5
    RPF;
    PoFBasedInspectionIntervalCalWithInspectionValues;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon1ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon1ActualCrack(k-1)];
end
RPF;
save coupon1_20221002_1.mat
toc;

%% 试件2
clear;clc;load('firstInspectionTest20221001.mat');
%
tic
inspectionValues=coupon2ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3;
trueValues=coupon2ActualCrack(k-1);
detectionNumber=10;

while inspectionValues(end)<1.27
    PFWithoutInspectionValues;
    PoFBasedInspectionIntervalCal;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon2ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon2ActualCrack(k-1)];
end
while inspectionValues(end)<2
    RPF;
    PoFBasedInspectionIntervalCalWithInspectionValues;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon2ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon2ActualCrack(k-1)];
end
RPF;
save coupon2_20221002_2.mat
toc;

%
% Intelligent inspection
% 试件3

clear;clc;load('firstInspectionTest20221001.mat');
%
tic
inspectionValues=coupon3ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3;
trueValues=coupon3ActualCrack(k-1);
detectionNumber=10;

while inspectionValues(end)<1.27
    PFWithoutInspectionValues;
    PoFBasedInspectionIntervalCal;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon3ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon3ActualCrack(k-1)];
end
while inspectionValues(end)<2
    RPF;
    PoFBasedInspectionIntervalCalWithInspectionValues;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon3ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon3ActualCrack(k-1)];
end
RPF;
save coupon3_20221002_2.mat
toc;


%%
% Intelligent inspection
% 试件4
clear;clc;load('firstInspectionTest20220921.mat');
%
tic
inspectionValues=coupon4ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3;
trueValues=coupon4ActualCrack(k-1);
detectionNumber=10;

while inspectionValues(end)<1.27
    PFWithoutInspectionValues;
    PoFBasedInspectionIntervalCal;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon4ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon4ActualCrack(k-1)];
end
while inspectionValues(end)<2.5
    RPF;
    PoFBasedInspectionIntervalCalWithInspectionValues;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon4ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon4ActualCrack(k-1)];
end
RPF;
save coupon4_20220925_1.mat
toc;


%% 试件6
% Intelligent inspection
% 试件6

clear;clc;load('firstInspectionTest20221001.mat');
%
tic
inspectionValues=coupon6ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3;
trueValues=coupon6ActualCrack(k-1);
detectionNumber=10;

while inspectionValues(end)<1.27
    PFWithoutInspectionValues;
    PoFBasedInspectionIntervalCal;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon6ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon6ActualCrack(k-1)];
end
%

while inspectionValues(end)<2
    RPF;
    PoFBasedInspectionIntervalCalWithInspectionValues;
    InspectionIntervalSet=[InspectionIntervalSet InspectionInterval];
    InspectionTime=sum(InspectionIntervalSet);
    probabilisticPropagation;
    inspectionValues=[inspectionValues coupon6ActualCrack(k-1)+abs(normrnd(0,0.1,1,1))/3];
    trueValues=[trueValues coupon6ActualCrack(k-1)];
end
RPF;
save coupon6_20221002_1.mat
toc;