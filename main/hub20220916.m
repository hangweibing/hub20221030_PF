%% ===================================================================
%% 程序初始化和环境设置
%% ===================================================================

% 清理工作空间
clear; close all; clc; tic; format long;

% 添加路径（可选）
% addpath surModelPackage

%% ===================================================================
%% 粒子滤波算法参数配置
%% ===================================================================

% 基本粒子滤波参数
n = 1;                                    % 状态向量的维度（每个粒子）
N = 1000;                                  % 粒子数量
v_sphere = 2;                             % 一维空间维度参数

% 正则化粒子滤波参数计算
A = (8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h = A*N^(-1/(n+4));                       % 正则化粒子滤波的平滑参数

%% ===================================================================
%% 数据存储矩阵初始化
%% ===================================================================

% 粒子滤波核心矩阵
D = zeros(44, 1000);                     % 协方差矩阵的对角线元素
e = zeros(N, 44, 1000);                  % Epanechikov核函数中的扰动项
cyclesperhour = 1950.70866;              % 每小时循环次数
Xpf = zeros(44, 1000);                   % 每一时刻的滤波估计值
xparticle = zeros(N, 44, 1000);          % 粒子状态矩阵
xparticle1 = zeros(N, 44, 1000);         % 重采样后的粒子
xparticle_cov = zeros(44, 44, 1000);     % 粒子协方差矩阵

% 参数粒子存储
upcrackparticles = zeros(N, 1000);       % 上表面裂纹粒子
logCstarparticles = zeros(N, 1000);      % logC*参数粒子
gammaparticles = zeros(N, 1000);         % gamma参数粒子
weight = zeros(N, 1000);                 % 粒子权重

%% ===================================================================
%% 观测数据配置
%% ===================================================================

% 观测时间点（单位：小时）
% t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
% z      =[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;

% 实际使用的观测数据（简化版）
t_check = [5.8741E+01 1.3869E+02  2.1810E+02 2.4204E+02 2.8120E+02];
z       = [2.4882E+00 4.0190E+00  1.3130E+01 1.6057E+01 2.0226E+01] + 30;

% 观测相关参数
zPred = zeros(N, 1);                     % 预测观测值
R = 0.5;                                 % 观测噪声方差

%% ===================================================================
%% 裂纹几何模型参数
%% ===================================================================

% 几何离散化参数
n_nodes = 21;                            % 节点数量
centers = [13, 30];                      % 圆心坐标
downCenter = [6, 37.82842712];           % 下圆圆心
downRadius = 3;                          % 下圆半径
nRegPoint = n_nodes;                     % 区域点数量

%% ===================================================================
%% 初始裂纹几何形状计算
%% ===================================================================

% 角度离散化（1/4圆弧：从270°到360°）
thetas = linspace(3/2*pi, 2*pi, n_nodes); % 角度分布

% 初始化坐标存储数组
yIniRegSet = zeros(1, n_nodes);          % 初始y坐标集
zIniRegSet = zeros(1, n_nodes);          % 初始z坐标集

%% ===================================================================
%% Walker定律参数先验分布设置
%% ===================================================================

% 先验分布参数（已注释的多元正态分布）
% mu = [-10.8 3];
% SIGMA = [0.05^2 -0.02^2; -0.02^2 0.01^2];

%% ===================================================================
%% 粒子群初始化
%% ===================================================================

for i = 1:N   % 遍历所有粒子
    %% 参数采样（Walker定律参数）
    % 生成参数样本（先验分布）- 已注释
    % paraSet = mvnrnd(mu, SIGMA, 1);
    % logCstar = paraSet(1);
    % gamma = paraSet(2);

    % 实际使用的采样策略
    logCstar = unifrnd(-11.2, -10.6, 1, 1);  % 均匀分布采样
    gamma = normrnd(3, 0.03, 1, 1);          % 正态分布采样

    %% 几何参数采样（裂纹尺寸）
    a = normrnd(2, 0.05, 1, 1);  % 裂纹半径（均值2mm，标准差0.05mm）
    c = a;                        % 设置为圆形（c = a）

    %% 裂纹轮廓坐标计算
    for j = 1:n_nodes
        % 参数方程：圆心 + 半径 × 方向向量
        yIniRegSet(j) = centers(1) + c*sin(thetas(j));
        zIniRegSet(j) = centers(2) + a*cos(thetas(j));
    end

    %% 粒子状态存储
    % 完整状态向量：[y坐标集(21), z坐标集(21), logC*, gamma]
    xparticle(i, :, 1) = [yIniRegSet, zIniRegSet, logCstar, gamma];
end

%% ===================================================================
%% 初始时刻统计量计算
%% ===================================================================

% 提取各参数的粒子分布
upcrackparticles(:, 1) = xparticle(:, 42, 1);    % 上表面裂纹长度
logCstarparticles(:, 1) = xparticle(:, 43, 1);   % logC*参数
gammaparticles(:, 1) = xparticle(:, 44, 1);      % gamma参数

% 初始化粒子权重（均匀分布）
weight(:, 1) = 1/N * ones(N, 1);

% 计算初始滤波估计（各参数的均值）
Xpf(:, 1) = (mean(xparticle(:, :, 1)))';

% 计算初始协方差矩阵
xparticle_cov(:, :, 1) = cov(xparticle(:, :, 1));

%% ===================================================================
%% 数据加载和预处理
%% ===================================================================
% 加载ASTM标准载荷谱数据
load('AsteixSpectraData.mat');

%% ===================================================================
%% POD模型数据加载
%% ===================================================================

load('pod_models.mat'); % 加载POD模型数据

%% ===================================================================
%% 频谱数据扩展处理
%% ===================================================================

% 扩展频谱数据以满足长时间计算需求（7倍扩展）
spectra = [spectra spectra spectra spectra spectra spectra spectra];

% 去除第一个数据点，确保数据对齐
spectra = spectra(2:end);

%% ===================================================================
%% 平均应力增量计算
%% ===================================================================

% 时间步长设置（每1000个应力循环作为一个时间步）
step = 1000;
aver_delta_sigma_set = zeros(1, ceil((length(spectra))/2/step));
k=1;

% 一个step内的平均 delta_sigma
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
%% ===================================================================
%% 模型参数配置
%% ===================================================================

% 裂纹阶段到模型编号的映射表
stage = [1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 8];

% 各阶段模型的测试误差集合（10个误差值）
testErrSet = [0.019803420755871 0.058377623733943 0.013343080193008 ...
              0.009221345729625 0.037283991994545 0.011408674471790 ...
              0.082711915490345 0.032375406062069 0.041104885753232 ...
              0.057544490601307];

%% ===================================================================
%% 粒子滤波主循环
%% ===================================================================

% 清屏显示
clc;

%% ===================================================================
%% 循环变量初始化
%% ===================================================================

% 分裂状态临时存储
SPLITTE_temp = zeros(1, N);

% 时间和观测计数器
m = 2;  % 时间步索引（从第2步开始，第1步为初始状态）
j = 1;  % 观测计数器

%% ===================================================================
%% 分裂模型数据预加载
%% ===================================================================

% 加载第3阶段分裂模型的数据
Uinput_splitted_1 = Uinput_splitted{1};
averInput_splitted_1 = averInput_splitted{1};

% 加载第5阶段分裂模型的数据
Uinput_splitted_2 = Uinput_splitted{2};
averInput_splitted_2 = averInput_splitted{2};

%% ===================================================================
%% 主时间推进循环
%% ===================================================================

while (m-1)*step/1950.70866 <= t_check(end)
    %% 时间步数据准备
    % 获取上一时刻的所有粒子状态
    xparticlem_1 = xparticle(:, :, m-1);

    % 获取当前时间段的平均应力增量
    aver_delta_sigma = aver_delta_sigma_set(m-1);

    %% ===================================================================
    %% 粒子预测步骤（对每个粒子进行状态更新）
    %% ===================================================================

    parfor (i = 1:N)
        %% 粒子级变量初始化
        curUinput = {};
        curAverInput = {};
        SPLITTED = SPLITTE_temp(i);  % 获取上一时刻的分裂状态
        
        %% 为每个并行线程设置独立的随机数流（parfor线程安全）
        % 原方案: 直接使用全局randi()，导致线程间竞态条件和不可重现结果
        % 问题: parfor中所有线程共享全局随机数生成器状态
        % 解决: 为每个线程创建独立随机流，避免线程间干扰
        % 使用粒子索引和时间步作为种子，确保可重现性
        stream = RandStream('mt19937ar', 'Seed', i + m*N);

        %% 粒子状态提取
        xparticlei = xparticlem_1(i, :);  % 当前粒子的完整状态向量

        %% 裂纹尺寸计算（用于阶段判断）
        a_up = xparticlei(42) - 30;     % 上表面裂纹长度（z坐标）
        a_down = xparticlei(22) - 30;   % 下表面裂纹长度（z坐标）

        %% 裂纹发展阶段判断
        % 220914: 从离散采样改为连续采样策略
        m_index = 0;

        while m_index == 0
            m_index = getModelIndexFunc(a_up, a_down); % 判断裂纹发展阶段

            if m_index == 0
                % 阶段判断失败：使用线程本地随机数重采样
                tmp_sel = randi(stream, [1, N], 1, 1);
                xparticlei = xparticlem_1(tmp_sel, :);
                a_up = xparticlei(42) - 30;
                a_down = xparticlei(22) - 30;
            end
        end

        %% ===================================================================
        %% 智能模型选择和数据配置
        %% ===================================================================

        if SPLITTED && (m_index==3 || m_index==5)
            %% 分裂模型分支（几何复杂情况）
            if m_index==3
                curUinput = Uinput_splitted_1;      % 第3阶段分裂POD基
                curAverInput = averInput_splitted_1; % 第3阶段分裂平均形状
            elseif m_index==5
                curUinput = Uinput_splitted_2;      % 第5阶段分裂POD基
                curAverInput = averInput_splitted_2; % 第5阶段分裂平均形状
            end
            m_name = sprintf('nn_stage%ds', stage(m_index));  % 分裂模型文件名
        else
            %% 整体模型分支（标准情况）
            curUinput = Uinput_integrated{stage(m_index)};       % 标准POD基
            curAverInput = averInput_integrated{stage(m_index)};  % 标准平均形状
            m_name = sprintf('nn_stage%d', stage(m_index));       % 整体模型文件名
        end

        %% ===================================================================
        %% 粒子状态分量提取
        %% ===================================================================

        % 几何坐标分量
        yRegSet = xparticlei(1:21);     % y坐标集（21个节点）
        zRegSet = xparticlei(22:42);    % z坐标集（21个节点）

        % Paris定律参数
        logCstar = xparticlei(43);       % 材料参数logC*
        gamma = xparticlei(44);          % 材料参数γ

        % 可选：绘制裂纹几何形状
        % figure; plot_geometry_20; plot(zRegSet, yRegSet); axis equal;

        %% 调用预测模型更新粒子状态
        [yRegSet, zRegSet, SPLITTED, logCstar, gamma] = ...
            a2aNew(yRegSet, zRegSet, aver_delta_sigma, m_name, ...
                   curUinput, curAverInput, logCstar, gamma, step, testErrSet);

        %% 更新粒子状态（parfor兼容：直接确保实数）
        % 原方案: xparticle(i, :, m) = [...]; xparticle = real(xparticle);
        % 问题: 第二行违反parfor切片规则，读取整个数组导致竞态条件
        % 解决: 在赋值时直接取实部，避免全局数组操作
        xparticle(i, :, m) = real([yRegSet, zRegSet, logCstar, gamma]);
        SPLITTE_temp(i) = SPLITTED;      % 更新分裂状态
    end

    %% 计算粒子滤波统计量
    weight(:, m) = 1/N * ones(N, 1);             % 均匀权重初始化
    Xpf(:, m) = (mean(xparticle(:, :, m)))';     % 状态均值估计
    xparticle_cov(:, :, m) = cov(xparticle(:, :, m)); % 状态协方差

    % 220913 已弃用，防止粒子群退化，但会导致去除过多粒子，因此保留所有粒子
    % tmp_outliner = isoutlier(xparticle(:,42,m));
    % if sum(tmp_outliner)~=0
    %     tmp_normal_particles=xparticle(~tmp_outliner,:,m);
    %     for k =1:N
    %         if tmp_outliner(k)==1
    %             tmp_sel=randi([1,size(tmp_normal_particles,1)],1,1);
    %             xparticle(k,:,m)=tmp_normal_particles(tmp_sel,:);
    %         end
    %     end
    % end

    %% 提取关键参数的历史记录
    upcrackparticles(:, m) = xparticle(:, 42, m);    % 上表面裂纹长度
    logCstarparticles(:, m) = xparticle(:, 43, m);   % logC*参数
    gammaparticles(:, m) = xparticle(:, 44, m);      % gamma参数
    % 可选：根据观测次数调整观测噪声
    % if j<=4; R=0.5; else R=0.3; end;

    %% 观测更新步骤（当到达观测时刻时）
    if ((m-1)*step/1950.70866 < t_check(j)) && (m*step/1950.70866 >= t_check(j))
        %% 计算似然权重
        for i = 1:N   % 观测更新
            zPred(i) = upcrackparticles(i, m);        % 预测观测值
            z1(i) = z(j) - zPred(i);                  % 观测残差
            weight(i, m) = inv(sqrt(2*pi*det(R))) * ...
                          exp(-0.5*(z1(i))*inv(R)*(z1(i))') + 1e-99;
        end

        %% 归一化权重
        weight(:, m) = weight(:, m) ./ sum(weight(:, m));

        %% 更新状态估计
        Xpf(:, m) = 0;
        for i = 1:N
            Xpf(:, m) = Xpf(:, m) + (weight(i, m) * xparticle(i, :, m))';
        end

        %% 更新协方差
        xparticle_cov(:, :, m) = 0;
        for i = 1:N
            xparticle_cov(:, :, m) = xparticle_cov(:, :, m) + ...
                weight(i, m) * (xparticle(i, :, m)' - Xpf(:, m)) * ...
                (xparticle(i, :, m)' - Xpf(:, m))';
        end

        %% 正则化重采样
        for i = 1:44
            D(i, m) = sqrt(xparticle_cov(i, i, m));      % 标准差
            e(:, i, m) = kernelsampling(N)';             % 核采样扰动
        end
        outindex = randomr(weight(:, m));                % 按权重重采样
        xparticle(:, :, m) = xparticle(outindex, :, m); 
        % 已注释：传统重采样方法
        % xparticle1(:,:,m)=xparticle(outindex,:,m);                                 % 重采样
        % for i=1:44
        %     xparticle(:,i,m)=xparticle1(:,i,m)+h*D(i,m)*e(:,i,m);              % 扰动
        %     xparticle(:,i,m)=rearrange(xparticle(:,i,m),xparticle1(:,i,m))';     % 参考文献：Dynamic Bayesian Network for Aircraft Wing Health Monitoring Digital Twin
        % end
        % for i=1:N
        %     [yNewSet,zNewSet,SPLITTE_temp(i)]=addConstraintNewSatgeFunc(xparticle(i,1:21,m),xparticle(i,22:42,m));
        %     [xparticle(i,1:21,m),xparticle(i,22:42,m),~] = crackRegular5Func(yNewSet,zNewSet,nRegPoint,'false');
        % end

        j = j + 1;  % 观测计数器递增
    end

    m = m + 1;  % 时间步递增
    disp(['已完成' num2str((m-1)*step/1950.70866) '小时，进行了' num2str(j-1) '次观测']);
end

%% ===================================================================
%% 结果后处理和可视化
%% ===================================================================
close all
figure

%% 数据处理
clear x y_0 y_1 y_2 y_3 y_33 y21 PoF3

% 计算统计量
for i = 1:m-1
    y_1(i) = prctile(upcrackparticles(:, i), 99.95) - 30;  % 99.95%分位数
    y_2(i) = prctile(upcrackparticles(:, i), 0.05) - 30;   % 0.05%分位数
    y_3(i) = Xpf(42, i) - 30;                              % 均值估计
    x(i) = (i-1) * step / 1950.70866;                      % 时间轴
end

%% 观测数据
t_check = [5.8741E+01 1.3869E+02  2.1810E+02 2.4204E+02 2.8120E+02];
z       = [2.4882E+00 4.0190E+00  1.3130E+01 1.6057E+01 2.0226E+01] + 30;

%% 绘制结果
plot(t_check, z-30, '^', 'linewidth', 5); hold on
plot(x, y_3, 'b', 'linewidth', 5); hold on;
plot(x, y_1, 'r--', 'linewidth', 5); hold on
plot(x, y_2, 'r--', 'linewidth', 5); hold off;


%% 图表格式设置
xlabel('Flight hours/h', 'FontSize', 30);
ylabel('Surface crack length/mm', 'FontSize', 30);
set(get(gca, 'xlabel'), 'fontname', 'Times New Roman');
set(get(gca, 'ylabel'), 'fontname', 'Times New Roman');
set(gca, 'fontname', 'Times New Roman');
set(gca, 'FontSize', 40);

% 图例设置
legend('   Experimental value', '   Prediction mean', '   99.9% bounds', 'FontSize', 40);
legend('boxoff')
legend('Location', 'best');

% 网格设置
grid on;
set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
% axis([0,300,1,21]);
% t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
% z=[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;
% plot(t_check,z-30,'v','linewidth',5);hold on
% 
% figure
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,1));plot(XI_ksdensity,F_ksdensity,'Color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,158));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,216));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,294));plot(XI_ksdensity,F_ksdensity,'b-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,353));plot(XI_ksdensity,F_ksdensity,'b-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,m-1));plot(XI_ksdensity,F_ksdensity,'r-.','LineWidth',5);hold on
% % [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,588));plot(XI_ksdensity,F_ksdensity,'r-','LineWidth',5);hold off
% xlabel('logC','FontSize',40);
% ylabel('PDF','FontSize',40);
% set(get(gca,'xlabel'),'fontweight','bold');
% set(get(gca,'ylabel'),'fontweight','bold');
% set(get(gca,'xlabel'),'fontname','Times New Roman');
% set(get(gca,'ylabel'),'fontname','Times New Roman');
% set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',30);
% legend('   Prior','   Inspection 1','   Inspection 2','   Inspection 3','   Inspection 4','   Inspection 5','   Inspection 6','FontSize',30);
% legend('Location','best');
% legend('boxoff')
% % axis([-10.94,-10.76,0,60]);
% 
% figure
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,1));plot(XI_ksdensity,F_ksdensity,'Color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,158));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,216));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,294));plot(XI_ksdensity,F_ksdensity,'b-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,353));plot(XI_ksdensity,F_ksdensity,'b-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,m-1));plot(XI_ksdensity,F_ksdensity,'r-.','LineWidth',5);hold on
% % [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,588));plot(XI_ksdensity,F_ksdensity,'r-','LineWidth',5);hold off
% xlabel('\gamma','FontSize',30);
% set(get(gca,'xlabel'),'fontweight','bold');
% ylabel('PDF','FontSize',30);
% set(get(gca,'xlabel'),'fontname','Times New Roman');
% set(get(gca,'ylabel'),'fontname','Times New Roman');
% set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',30);
% legend('   Prior','   Inspection 1','   Inspection 2','   Inspection 3','   Inspection 4','   Inspection 5','   Inspection 6','FontSize',30);
% legend('Location','best');
% legend('boxoff')
% % axis([2.85,3.12,0,190]);

%%
% % regularized particle filter
% for m=2:1850       %��ǰԤ��
%     xparticlem_1=xparticle(:,:,m-1);
%     aver_delta_sigma=aver_delta_sigma_set(m-1);
%     parfor i=1:N
%         curModel={};
%         curUinput={};
%         curAverInput={};
%         SPLITTED=SPLITTE_temp(i);
%         xparticlei=xparticlem_1(i,:);
%         a_up=xparticlei(42)-30;
%         a_down=xparticlei(22)-30;
%         m_index = getModelIndexFunc(a_up,a_down); % �ж����ƽ׶�
%         if ~SPLITTED % ���ǰԵδ���룬���������ƴ���ģ�ͼ�����ȡ����ģ��
%             curModel=model_integrated{m_index};
%             curUinput=Uinput_integrated{m_index};
%             curAverInput=averInput_integrated{m_index};
%         else % ���ǰԵ�ѷ��룬�ӷ������ƴ���ģ�ͼ�����ȡ����ģ��
%             if m_index==3
%                 curModel=model_splitted{1};
%                 curUinput=Uinput_splitted{1};
%                 curAverInput=averInput_splitted{1};
%             elseif m_index==5
%                 curModel=model_splitted{2};
%                 curUinput=Uinput_splitted{2};
%                 curAverInput=averInput_splitted{2};
%             end
%         end
%         yRegSet=xparticlei(1:21);
%         zRegSet=xparticlei(22:42);
%         logCstar=xparticlei(43);
%         gamma=xparticlei(44);
%         [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aFunc(yRegSet,zRegSet,aver_delta_sigma,curModel,curUinput,curAverInput,logCstar,gamma,step);
%         xparticle(i,:,m)=[yRegSet,zRegSet,logCstar,gamma];
%         SPLITTE_temp(i)=SPLITTED;
%     end
%     disp(m)
%     weight(:,m)=1/N*ones(N,1);
%     Xpf(:,m)=(mean(xparticle(:,:,m)))';
%     xparticle_cov(:,:,m)=cov(xparticle(:,:,m));
%     upcrackparticles(:,m)=xparticle(:,42,m);
%     logCstarparticles(:,m)=xparticle(:,43,m);
%     gammaparticles(:,m)=xparticle(:,44,m);
% end
% parfor i=1:N   %�������
%     zPred(i)=upcrackparticles(i,117);
%     z1(i)=z(1)-zPred(i);
%     weight(i,117)=inv(sqrt(2*pi*det(R)))*exp(-0.5*(z1(i))*inv(R)*(z1(i))')+1e-99;
% end
% weight(:,117)=weight(:,117)./sum(weight(:,117));
% Xpf(:,117)=0;
% for i=1:N
%     Xpf(:,117)=Xpf(:,117)+(weight(i,117)*xparticle(i,:,117))';
% end
% xparticle_cov(:,:,117)=0;
% for i=1:N
%     xparticle_cov(:,:,117)=xparticle_cov(:,:,117)+weight(i,117)*(xparticle(i,:,117)'-Xpf(:,117))*(xparticle(i,:,117)'-Xpf(:,117))';
% end
% for i=1:44
%     D(i,117)=sqrt(xparticle_cov(i,i,117));
%     e(:,i,117)=kernelsampling(N)';
% end
% outindex=randomr(weight(:,117));
% xparticle1(:,:,117)=xparticle(outindex,:,117);                                 %�ز���
% for i=1:44
%     xparticle(:,i,117)=xparticle1(:,i,117)+h*D(i,117)*e(:,i,117);              %����
%     xparticle(:,i,117)=rearrange(xparticle(:,i,117),xparticle1(:,i,117))';     %�μ����ס�Dynamic Bayesian Network for Aircraft Wing Health Monitoring Digital Twin��
% end
% for i=1:N
%     [yNewSet,zNewSet,SPLITTE_temp(i)]=addConstraintNewSatgeFunc(xparticle(i,1:21,117),xparticle(i,22:42,117));
%     [xparticle(i,1:21,117),xparticle(i,22:42,117),~] = crackRegular5Func(yNewSet,zNewSet,nRegPoint,'false');
% end
% 
% %�ڶ����غ���
% for m=2:81       %��ǰԤ��
%     xparticlem_1=xparticle(:,:,m+115);
%     aver_delta_sigma=aver_delta_sigma_set2(m-1);
%     M=m+116;
%     parfor i=1:N
%         curModel={};
%         curUinput={};
%         curAverInput={};
%         SPLITTED=SPLITTE_temp(i);
%         xparticlei=xparticlem_1(i,:);
%         a_up=xparticlei(42)-30;
%         a_down=xparticlei(22)-30;
%         m_index = getModelIndexFunc(a_up,a_down); % �ж����ƽ׶�
%         if ~SPLITTED % ���ǰԵδ���룬���������ƴ���ģ�ͼ�����ȡ����ģ��
%             curModel=model_integrated{m_index};
%             curUinput=Uinput_integrated{m_index};
%             curAverInput=averInput_integrated{m_index};
%         else % ���ǰԵ�ѷ��룬�ӷ������ƴ���ģ�ͼ�����ȡ����ģ��
%             if m_index==3
%                 curModel=model_splitted{1};
%                 curUinput=Uinput_splitted{1};
%                 curAverInput=averInput_splitted{1};
%             elseif m_index==5
%                 curModel=model_splitted{2};
%                 curUinput=Uinput_splitted{2};
%                 curAverInput=averInput_splitted{2};
%             end
%         end
%         yRegSet=xparticlei(1:21);
%         zRegSet=xparticlei(22:42);
%         logCstar=xparticlei(43);
%         gamma=xparticlei(44);
%         [yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aFunc(yRegSet,zRegSet,aver_delta_sigma,curModel,curUinput,curAverInput,logCstar,gamma,step);
%         xparticle(i,:,M)=[yRegSet,zRegSet,logCstar,gamma];
%         SPLITTE_temp(i)=SPLITTED;
%     end
%     m/81
%     2
%     weight(:,M)=1/N*ones(N,1);
%     Xpf(:,M)=(mean(xparticle(:,:,M)))';
%     xparticle_cov(:,:,M)=cov(xparticle(:,:,M));
%     upcrackparticles(:,M)=xparticle(:,42,M);
%     logCstarparticles(:,M)=xparticle(:,43,M);
%     gammaparticles(:,M)=xparticle(:,44,M);
% end
% parfor i=1:N   %�������
%     zPred(i)=upcrackparticles(i,197);
%     z1(i)=z(2)-zPred(i);
%     weight(i,197)=inv(sqrt(2*pi*det(R)))*exp(-0.5*(z1(i))*inv(R)*(z1(i))')+1e-99;
% end
% weight(:,197)=weight(:,197)./sum(weight(:,197));
% Xpf(:,197)=0;
% for i=1:N
%     Xpf(:,197)=Xpf(:,197)+(weight(i,197)*xparticle(i,:,197))';
% end
% xparticle_cov(:,:,197)=0;
% for i=1:N
%     xparticle_cov(:,:,197)=xparticle_cov(:,:,197)+weight(i,197)*(xparticle(i,:,197)'-Xpf(:,197))*(xparticle(i,:,197)'-Xpf(:,197))';
% end
% for i=1:44
%     D(i,197)=sqrt(xparticle_cov(i,i,197));
%     e(:,i,197)=kernelsampling(N)';
% end
% outindex=randomr(weight(:,197));
% xparticle1(:,:,197)=xparticle(outindex,:,197);                             %�ز���
% for i=1:44
%     xparticle(:,i,197)=xparticle1(:,i,197)+h*D(i,197)*e(:,i,197);          %����
%     xparticle(:,i,197)=rearrange(xparticle(:,i,197),xparticle1(:,i,197))'; %�μ����ס�Dynamic Bayesian Network for Aircraft Wing Health Monitoring Digital Twin��
% end
% for i=1:N
%     [yNewSet,zNewSet,SPLITTE_temp(i)]=addConstraintNewSatgeFunc(xparticle(i,1:21,197),xparticle(i,22:42,197));
%     [xparticle(i,1:21,197),xparticle(i,22:42,197),~] = crackRegular5Func(yNewSet,zNewSet,nRegPoint,'false');
% end
% 
% 
% % %% post processing
% % close all;clc;clear all;load('hub20211123.mat')
% %
% % logCstar_post_Estimation=[Xpf(43,117) Xpf(43,197) Xpf(43,274) Xpf(43,338) Xpf(43,360)];
% % gamma_post_Estimation=[Xpf(44,117) Xpf(44,197) Xpf(44,274) Xpf(44,338) Xpf(44,360)];
% clc;figure
% % [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,1));plot(XI_ksdensity,F_ksdensity,'LineWidth',2);hold on
% plot([sort(xparticle(:,43,1));-10.8],[0 ones(1,(length(xparticle(:,43,1))-1)).*10 0],'Color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,117));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,197));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,274));plot(XI_ksdensity,F_ksdensity,'b-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,338));plot(XI_ksdensity,F_ksdensity,'b-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,43,360));plot(XI_ksdensity,F_ksdensity,'r-','LineWidth',5);hold off
% xlabel('logC','FontSize',40);
% ylabel('PDF','FontSize',40);
% % set(get(gca,'xlabel'),'fontweight','bold');
% % set(get(gca,'ylabel'),'fontweight','bold');
% set(get(gca,'xlabel'),'fontname','Times New Roman');
% set(get(gca,'ylabel'),'fontname','Times New Roman');
% set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',30);
% legend('   Prior','   t = 59.39h','   t = 100.34h','   t = 139.59h','   t = 172.35h','   t = 183.62h','FontSize',30);
% legend('Location','best');
% legend('boxoff')
% axis([-10.94,-10.76,0,60]);
% %2.4882E+00
% 
% figure
% plot(sort(xparticle(:,44,1)),[0 ones(1,(length(xparticle(:,44,1))-2)).*5 0],'Color',[0.5 0.5 0.5],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,117));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,197));plot(XI_ksdensity,F_ksdensity,'Color',[0 0 0],'LineStyle','-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,274));plot(XI_ksdensity,F_ksdensity,'b-.','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,338));plot(XI_ksdensity,F_ksdensity,'b-','LineWidth',5);hold on
% [F_ksdensity,XI_ksdensity]=ksdensity(xparticle(:,44,360));plot(XI_ksdensity,F_ksdensity,'r-','LineWidth',5);hold off
% xlabel('\gamma','FontSize',30);
% set(get(gca,'xlabel'),'fontweight','bold');
% ylabel('PDF','FontSize',30);
% set(get(gca,'xlabel'),'fontname','Times New Roman');
% set(get(gca,'ylabel'),'fontname','Times New Roman');
% set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',30);
% legend('   Prior','   t = 59.39h','   t = 100.34h','   t = 139.59h','   t = 172.35h','   t = 183.62h','FontSize',30);
% legend('Location','best');
% legend('boxoff')
% axis([2.85,3.12,0,190]);
% % figure(3)
% % hist(RUL,10)
% % legend('RUL2','FontSize',20);
% %
% figure
% timeSeries_1=zeros(1,ceil((length(spectra1))/2/step));
% for i=1:floor((length(spectra1))/2/step)
%     timeSeries_1(i)=i*step/1950.70866;
% end
% timeSeries_1(end)=t_check(1);
% 
% % timeSeries_2=zeros(1,ceil((length(spectra2))/2/step));
% % for i=1:floor((length(spectra2))/2/step)
% %     timeSeries_2(i)=t_check(1)+i*step/1950.70866;
% % end
% % timeSeries_2(end)=t_check(2);
% % 
% % timeSeries_3=zeros(1,ceil((length(spectra3))/2/step));
% % for i=1:floor((length(spectra3))/2/step)
% %     timeSeries_3(i)=t_check(2)+i*step/1950.70866;
% % end
% % timeSeries_3(end)=t_check(3);
% % 
% % timeSeries_4=zeros(1,ceil((length(spectra4))/2/step));
% % for i=1:floor((length(spectra4))/2/step)
% %     timeSeries_4(i)=t_check(3)+i*step/1950.70866;
% % end
% % timeSeries_4(end)=t_check(4);
% % 
% % timeSeries_5=zeros(1,ceil((length(spectra5))/2/step));
% % for i=1:floor((length(spectra5))/2/
% t=[0 timeSeries_1];
% 
% clear x y y1 y2
% for i=1:321
%     x(i)=t(i);
%     %     y(i)=prctile(upcrackparticles(:,i),50);
%     y(i)=Xpf(42,i)-30;
% end
% plot(x,y,'b','linewidth',5);hold on
% for i=1:321
%     x(i)=t(i);
%     y1(i)=prctile(upcrackparticles(:,i),97.5)-30;
% end
% plot(x,y1,'r--','linewidth',5);hold on
% for i=1:321
%     x(i)=t(i);
%     y2(i)=prctile(upcrackparticles(:,i),2.5)-30;
% end
% plot(x,y2,'r--','linewidth',5);
% 
% xlabel('Flight hours/h','FontSize',30);
% ylabel('Surface crack length/mm','FontSize',30);
% set(get(gca,'xlabel'),'fontname','Times New Roman');
% set(get(gca,'ylabel'),'fontname','Times New Roman');
% set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',30);
% legend('   Prediction mean', '   95% bounds','FontSize',30);
% legend('boxoff')
% legend('Location','North');
% grid on;
% set(gca,'gridlinestyle',':','gridcolor','k');