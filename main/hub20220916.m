%% ===================================================================
%% 程序初始化和环境设置
%% ===================================================================

% 清理工作空间
clear; close all; clc; tic; format long;

% 设置图形参数
set(groot, 'defaultFigureColor', 'white');
set(groot, 'defaultAxesFontSize', 10);
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');

%% ===================================================================
%% 随机数种子设置（确保结果可重复）
%% ===================================================================
SIM_SEED = 2023;  % 您可以修改此数值来获得不同的随机序列
rng(SIM_SEED);    % 设置全局随机数种子

% 添加路径（可选）
% addpath surModelPackage

% 初始化 KDE 采样
kde_sampling('preprocess');


%% ===================================================================
%% 粒子滤波算法参数配置
%% ===================================================================

% Debug模式开关
DEBUG_MODE = false;                    % 设置为 true 开启debug模式，false 关闭

% 基本粒子滤波参数
n = 1;                                    % 状态向量的维度（每个粒子）
N = 10000;                                  % 粒子数量
v_sphere = 2;                             % 一维空间维度参数

% 正则化粒子滤波参数计算
A = (8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h = A*N^(-1/(n+4));                       % 正则化粒子滤波的平滑参数

% 预采样所有粒子的参数组合（提高效率）
if ~DEBUG_MODE
    fprintf('预采样 %d 个粒子的参数组合...\n', N);
    particle_params = kde_sampling('batch_sample', N);
    % particle_params = copula_sampling('batch_sample', N);
    fprintf('参数预采样完成！\n');
else
    fprintf('[Debug模式] 跳过预采样，将使用训练数据直接采样\n');
    particle_params = [];  % Debug模式下不需要预采样
end

%% ===================================================================
%% 采样结果可视化
%% ===================================================================
if ~DEBUG_MODE
    % 获取训练数据用于可视化对比
    T = readtable('AM-TC4-GRO.xlsx');
    X_train = table2array(T);
    param_names = T.Properties.VariableNames;
    d = size(X_train, 2);

    % 创建简单的可视化对比
    figure('Name', '采样参数分布对比', 'Position', [100, 100, 1000, 600]);
    for i = 1:d
        subplot(2,2,i)
        % 绘制采样数据的分布
        histogram(particle_params(:,i), 30, 'Normalization','pdf', 'FaceColor', 'b', 'FaceAlpha', 0.6);
        hold on
        % 绘制训练数据的分布
        histogram(X_train(:,i), 15, 'Normalization','pdf', 'FaceColor', 'r', 'FaceAlpha', 0.4);

        % 设置标题和标签
        clean_name = strrep(param_names{i}, '_', '-');
        title(sprintf('%s Distribution Comparison', clean_name), 'Interpreter', 'none', 'FontSize', 11, 'FontName', 'Arial')
        xlabel(clean_name, 'Interpreter', 'none', 'FontSize', 9, 'FontName', 'Arial')
        ylabel('Probability Density', 'FontSize', 9, 'FontName', 'Arial')
        legend(' Sampled', 'Training Data', 'Location', 'best', 'FontSize', 8, 'FontName', 'Arial')
        grid on
        set(gca, 'FontName', 'Arial', 'FontSize', 8)
    end

    %% 三参数联合分布可视化（三维散点图）
    figure('Name', '三参数联合分布对比', 'Position', [100, 100, 1200, 800]);

    % 计算需要显示的三参数组合 (4个参数的三参数组合 = 4种)
    param_triplets = [
        1,2,3;  % 参数1,2,3
        1,2,4;  % 参数1,2,4
        1,3,4;  % 参数1,3,4
        2,3,4   % 参数2,3,4
        ];

    for k = 1:4
        i = param_triplets(k,1);
        j = param_triplets(k,2);
        l = param_triplets(k,3);

        subplot(2,2,k);

        % 绘制训练数据的三维散点图
        scatter3(X_train(:,i), X_train(:,j), X_train(:,l), 20, 'filled', ...
            'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
        hold on;

        % 绘制Copula采样数据的三维散点图
        scatter3(particle_params(:,i), particle_params(:,j), particle_params(:,l), 15, 'filled', ...
            'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');

        % 清理变量名用于显示（移除下划线等特殊字符）
        clean_name_i = strrep(param_names{i}, '_', '-');
        clean_name_j = strrep(param_names{j}, '_', '-');
        clean_name_l = strrep(param_names{l}, '_', '-');

        xlabel(clean_name_i, 'Interpreter', 'none', 'FontSize', 9, 'FontName', 'Arial');
        ylabel(clean_name_j, 'Interpreter', 'none', 'FontSize', 9, 'FontName', 'Arial');
        zlabel(clean_name_l, 'Interpreter', 'none', 'FontSize', 9, 'FontName', 'Arial');

        title(sprintf('%s vs %s vs %s', clean_name_i, clean_name_j, clean_name_l), ...
            'Interpreter', 'none', 'FontSize', 10, 'FontName', 'Arial');

        legend('Training Data', 'Sampled', 'Location', 'best', 'FontSize', 8, 'FontName', 'Arial');

        grid on;
        set(gca, 'FontName', 'Arial', 'FontSize', 8);

        % 设置视角以获得更好的视觉效果
        view(45, 30);
    end

    fprintf('参数分布可视化完成！\n');
end
%% ===================================================================
%% 数据存储矩阵初始化
%% ===================================================================

% 粒子滤波核心矩阵
D = zeros(46, 1000);                     % 协方差矩阵的对角线元素（扩展为46维）
e = zeros(N, 46, 1000);                  % Epanechikov核函数中的扰动项（扩展为46维）
cyclesperhour = 1950.70866;              % 每小时循环次数
Xpf = zeros(46, 1000);                   % 每一时刻的滤波估计值（扩展为46维）
xparticle = zeros(N, 46, 1000);          % 粒子状态矩阵（扩展为46维以容纳4个NASGRO参数）
xparticle1 = zeros(N, 46, 1000);         % 重采样后的粒子
xparticle_cov = zeros(46, 46, 1000);     % 粒子协方差矩阵

% 参数粒子存储
upcrackparticles = zeros(N, 1000);       % 上表面裂纹粒子
log_theta1_particles = zeros(N, 1000);   % log_theta1_参数粒子（NASGRO模型）
theta2_particles = zeros(N, 1000);       % theta2参数粒子（NASGRO模型）
theta3_particles = zeros(N, 1000);       % theta3参数粒子（NASGRO模型）
k2_particles = zeros(N, 1000);           % k2参数粒子（NASGRO模型）
weight = zeros(N, 1000);                 % 粒子权重

% POF计算相关变量
particles_K_max = zeros(N, 1);          % 每个粒子的最大应力强度因子
POF_array = zeros(1000, 1);              % 失效概率数组（每个时间步）
mu_Kc = 33.4 ;                            % 断裂韧性均值 (MPa√m)
std_Kc = 3.34 ;                           % 断裂韧性标准差 (MPa√m)

%% ===================================================================
%% 观测数据配置
%% ===================================================================

% 观测时间点（单位：小时）
% t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
% z      =[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;

% 实际使用的观测数据（铝合金）
% t_check = [9.9534E+01 1.3869E+02  1.8275E+02 2.4204E+02 2.8120E+02];
% z       = [2.9894E+00 4.0190E+00  7.2212E+00 1.6057E+01 2.0226E+01] + 30;

% 实际使用的观测数据（钛合金）
t_check = [1.6869E+02 2.1275E+02  2.7204E+02 3.1120E+02 3.4546E+02];
z       = [2.9894E+00 4.0190E+00  7.2212E+00 1.6057E+01 2.0226E+01] + 30;

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
%% 粒子群初始化
%% ===================================================================

for i = 1:N   % 遍历所有粒子
    %% 参数采样
    if DEBUG_MODE
        %% Debug模式：使用训练数据直接采样
        [log_theta1_, theta2, theta3, k2] = debug_sampling();
    else
        %% 正常模式：使用预采样的 Copula 参数
        log_theta1_ = particle_params(i, 1);
        theta2 = particle_params(i, 2);
        theta3 = particle_params(i, 3);
        k2 = particle_params(i, 4);
    end

    % 使用新的变量名（已从预采样参数映射）
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
    % 完整状态向量：[y坐标集(21), z坐标集(21), log_theta1_, theta2, theta3, k2]
    xparticle(i, :, 1) = [yIniRegSet, zIniRegSet, log_theta1_, theta2, theta3, k2];
end

%% ===================================================================
%% 初始时刻统计量计算
%% ===================================================================

% 提取各参数的粒子分布
upcrackparticles(:, 1) = xparticle(:, 42, 1);      % 上表面裂纹长度
log_theta1_particles(:, 1) = xparticle(:, 43, 1);  % log_theta1_参数（NASGRO模型）
theta2_particles(:, 1) = xparticle(:, 44, 1);      % theta2参数（NASGRO模型）
theta3_particles(:, 1) = xparticle(:, 45, 1);      % theta3参数（NASGRO模型）
k2_particles(:, 1) = xparticle(:, 46, 1);          % k2参数（NASGRO模型）

% 初始化粒子权重（均匀分布）
weight(:, 1) = 1/N * ones(N, 1);

% 初始化粒子坐标信息存储（cell数组，1×N，每列保存一个粒子的坐标历史）
particles_coordinates = cell(1, N);

% 计算初始滤波估计（各参数的均值）
Xpf(:, 1) = (mean(xparticle(:, :, 1)))';

% 计算初始协方差矩阵
xparticle_cov(:, :, 1) = cov(xparticle(:, :, 1));

%% ===================================================================
%% 数据加载和预处理
%% ===================================================================
% 加载ASTM标准载荷谱数据
load('AsteixSpectraData.mat');
load('pod_models.mat'); % 加载POD模型数据

spectra = [spectra spectra spectra spectra spectra spectra spectra];
% 去除第一个数据点
spectra = spectra(2:end);

%% ===================================================================
%% 平均应力增量计算
%% ===================================================================

% 时间步长设置（每1000个应力循环作为一个时间步）
step = 1000;
aver_delta_sigma_set = zeros(1, ceil((length(spectra))/2/step));
aver_R_set = zeros(1, ceil((length(spectra))/2/step));  % 平均应力比数组
aver_Smax_set = zeros(1, ceil((length(spectra))/2/step)); % 平均最大应力数组
k=1;

% 一个step内的平均 delta_sigma 和平均应力比
for i=1:floor((length(spectra))/2/step)
    delta_sigmas=zeros(1,step);
    Smax_values=zeros(1,step);
    Smin_values=zeros(1,step);
    for j=1:step
        Smax=spectra(2*(k+j-1));
        Smin=spectra(2*(k+j-1)-1);
        delta_sigmas(j)=Smax-Smin;
        Smax_values(j)=Smax;
        Smin_values(j)=Smin;
    end
    aver_delta_sigma_set(i)=mean(delta_sigmas);
    % 计算平均应力比 R = Smin的均值 / Smax的均值
    aver_R_set(i) = mean(Smin_values) / mean(Smax_values);
    % 计算平均最大应力
    aver_Smax_set(i) = mean(Smax_values);
    k=k+step;
end
%


delta_sigmas=0;
Smax_values=0;
Smin_values=0;
k=floor((length(spectra))/2/step)*step+1;
i=1;
while k<=(length(spectra)/2)
    Smax=spectra(2*k);
    Smin=spectra(2*k-1);
    delta_sigmas(i)=Smax-Smin;
    Smax_values(i)=Smax;
    Smin_values(i)=Smin;
    k=k+1;
    i=i+1;
end
aver_delta_sigma_set(end)=mean(delta_sigmas);
aver_R_set(end) = mean(Smin_values) / mean(Smax_values);
aver_Smax_set(end) = mean(Smax_values);

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
%% Debug模式配置（裂纹可视化）
%% ===================================================================

% Debug模式参数（仅在DEBUG_MODE=true时生效）
DEBUG_PARTICLE_IDX = 228;               % 需要跟踪绘制的粒子编号（1~N）
DEBUG_PLOT_INTERVAL = 10;             % 绘图间隔：每隔多少次while循环保存一次图像

% Debug模式计数器初始化
if DEBUG_MODE
    debug_plot_counter = 0;           % 绘图计数器
    fprintf('[Debug模式] 已开启，将跟踪粒子 %d 的裂纹变化，绘图间隔：每 %d 次循环\n', ...
        DEBUG_PARTICLE_IDX, DEBUG_PLOT_INTERVAL);
end

%% ===================================================================
%% 粒子滤波主循环
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

total_tic = tic;  % 初始化总耗时计时器
while (m-1)*step/1950.70866 <= t_check(end)
    iter_tic = tic;  % 初始化当前步耗时计时器
    %% 时间步数据准备
    % 获取上一时刻的所有粒子状态
    xparticlem_1 = xparticle(:, :, m-1);

    % 获取当前时间段的平均应力增量
    aver_delta_sigma = aver_delta_sigma_set(m-1);

    %% ===================================================================
    %% 粒子预测步骤（对每个粒子进行状态更新）
    %% ===================================================================

    % 初始化当前时间步的K值存储
    particles_deltaK_max_temp = zeros(N, 1);

    parfor (i = 1:N)
        try
            %% 粒子级变量初始化
            curUinput = {};
            curAverInput = {};
            SPLITTED = SPLITTE_temp(i);  % 获取上一时刻的分裂状态

            %% 为每个并行线程设置独立的随机数流（parfor线程安全）
            stream = RandStream('mt19937ar', 'Seed', SIM_SEED + i + m*N);

            %% 粒子状态提取
            xparticlei = xparticlem_1(i, :);  % 当前粒子的完整状态向量

            % 如果上一时刻已经是NaN，则直接跳过
            if any(isnan(xparticlei))
                error('A2A:InvalidParticle', 'Inherited NaN');
            end

            %% 裂纹尺寸计算（用于阶段判断）
            a_up = xparticlei(42) - 30;     % 上表面裂纹长度（z坐标）
            a_down = xparticlei(22) - 30;   % 下表面裂纹长度（z坐标）

            %% 裂纹发展阶段判断
            m_index = 0;
            while m_index == 0
                m_index = getModelIndexFunc(a_up, a_down);

                if m_index == 0
                    tmp_sel = randi(stream, [1, N], 1, 1);
                    xparticlei = xparticlem_1(tmp_sel, :);
                    if any(isnan(xparticlei)), error('A2A:InvalidParticle', 'Selected NaN'); end
                    a_up = xparticlei(42) - 30;
                    a_down = xparticlei(22) - 30;
                end
            end

            %% 模型选择
            if SPLITTED && (m_index==3 || m_index==5)
                if m_index==3
                    curUinput = Uinput_splitted_1;
                    curAverInput = averInput_splitted_1;
                elseif m_index==5
                    curUinput = Uinput_splitted_2;
                    curAverInput = averInput_splitted_2;
                end
                m_name = sprintf('nn_stage%ds', stage(m_index));
            else
                curUinput = Uinput_integrated{stage(m_index)};
                curAverInput = averInput_integrated{stage(m_index)};
                m_name = sprintf('nn_stage%d', stage(m_index));
            end

            %% 粒子状态分量提取
            yRegSet = xparticlei(1:21);
            zRegSet = xparticlei(22:42);
            particles_coordinates{i} = [particles_coordinates{i}; [yRegSet, zRegSet]];

            log_theta1_ = xparticlei(43);
            theta2 = xparticlei(44);
            theta3 = xparticlei(45);
            k2 = xparticlei(46);

            %% 调用预测模型更新粒子状态
            [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, deltaKSet] = ...
                a2aNew(yRegSet, zRegSet, aver_delta_sigma, aver_R_set(m-1), aver_Smax_set(m-1), m_name, ...
                curUinput, curAverInput, log_theta1_, theta2, theta3, k2, step, testErrSet, i);

            particles_deltaK_max_temp(i) = deltaKSet(end);
            xparticle(i, :, m) = real([yRegSet, zRegSet, log_theta1_, theta2, theta3, k2]);
            SPLITTE_temp(i) = SPLITTED;

        catch ME
            if strcmp(ME.identifier, 'A2A:InvalidParticle')
                % 捕获到无效粒子，标记为 NaN 并不参与后续统计
                xparticle(i, :, m) = NaN;
                particles_deltaK_max_temp(i) = NaN;
                SPLITTE_temp(i) = 0;
            else
                rethrow(ME);
            end
        end
    end

    % Debug模式：统一保存所有粒子的坐标信息
    if DEBUG_MODE
        save('debug_particles_coordinates.mat', 'particles_coordinates', 'm');
    end

    % 计算粒子滤波统计量 (排除 NaN 无效粒子)
    valid_mask = ~isnan(xparticle(:, 1, m));
    weight(:, m) = 1/N * ones(N, 1);
    if any(valid_mask)
        Xpf(:, m) = (mean(xparticle(valid_mask, :, m)))';
        xparticle_cov(:, :, m) = cov(xparticle(valid_mask, :, m));
    else
        Xpf(:, m) = Xpf(:, m-1); % 如果全部失效，保持不变
    end

    %% ===================================================================
    %% 计算POF（！！！！！注意K的单位换算！！！！！）
    %% ===================================================================

    try
        % deltaK_max 计算出 K_max
        valid_K = particles_deltaK_max_temp(~isnan(particles_deltaK_max_temp));
        if ~isempty(valid_K)
            particles_K_max_valid = valid_K / (1 - aver_R_set(m-1));
            POF_array(m) = calculatePOF(particles_K_max_valid / sqrt(1000), mu_Kc, std_Kc);
        else
            POF_array(m) = POF_array(m-1);
        end
    catch ME
        warning('POF计算失败 (时间步 %d): %s', m, ME.message);
        POF_array(m) = 0;  % 失败时设为0
    end

    %% 提取关键参数的历史记录
    upcrackparticles(:, m) = xparticle(:, 42, m);      % 上表面裂纹长度
    log_theta1_particles(:, m) = xparticle(:, 43, m);  % log_theta1_参数（NASGRO模型）
    theta2_particles(:, m) = xparticle(:, 44, m);      % theta2参数（NASGRO模型）
    theta3_particles(:, m) = xparticle(:, 45, m);      % theta3参数（NASGRO模型）
    k2_particles(:, m) = xparticle(:, 46, m);          % k2参数（NASGRO模型）
    % 可选：根据观测次数调整观测噪声
    % if j<=4; R=0.5; else R=0.3; end;

    %% 观测更新步骤（当到达观测时刻时）
    if ((m-1)*step/1950.70866 < t_check(j)) && (m*step/1950.70866 >= t_check(j))
        %% 计算似然权重
        for i = 1:N   % 观测更新
            zPred(i) = upcrackparticles(i, m);
            if isnan(zPred(i))
                weight(i, m) = 1e-99; % 对无效粒子赋予极小权重
                continue;
            end
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
        % outindex = randomr(weight(:, m));                % 按权重重采样
        outindex = randomr(weight(:, m));  % 按权重重采样，使用确定性种子
        xparticle(:, :, m) = xparticle(outindex, :, m);

        % 重采样后清空粒子坐标历史，重新开始记录
        particles_coordinates = cell(1, N);

        j = j + 1;  % 观测计数器递增
    end

    %% ===================================================================
    %% Debug模式：裂纹可视化绘图
    %% ===================================================================
    if DEBUG_MODE
        debug_plot_counter = debug_plot_counter + 1;

        % 每隔指定次数循环，绘制指定粒子的裂纹图像
        if debug_plot_counter >= DEBUG_PLOT_INTERVAL
            % 获取指定粒子的当前裂纹坐标
            yDebugSet = xparticle(DEBUG_PARTICLE_IDX, 1:21, m);
            zDebugSet = xparticle(DEBUG_PARTICLE_IDX, 22:42, m);

            % 计算当前飞行小时数
            current_flight_hours = (m-1)*step/1950.70866;

            % 调用绘图函数
            plotCrackCoordinates(yDebugSet, zDebugSet, DEBUG_PARTICLE_IDX, m, current_flight_hours);
            fprintf('  [Debug] 已绘制粒子 %d 在时间步 %d (%.2f小时) 的裂纹图像\n', ...
                DEBUG_PARTICLE_IDX, m, current_flight_hours);


            % 重置计数器
            debug_plot_counter = 0;
        end
    end

    m = m + 1;  % 时间步递增

    % 计算耗时
    iter_time = toc(iter_tic);
    total_time = toc(total_tic);

    disp(['已完成' num2str((m-1)*step/1950.70866) '小时，进行了' num2str(j-1) '次观测', ...
        '，当前步耗时：' num2str(iter_time, '%.2f') 's', ...
        '，累计总耗时：' num2str(total_time, '%.2f') 's']);
end

%% ===================================================================
%% 结果后处理和可视化
%% ===================================================================
close all
%% 数据处理
clear x y_0 y_1 y_2 y_3 y_33 y21 PoF3

% 计算统计量
for i = 1:m-1
    valid_up = upcrackparticles(:, i);
    valid_up = valid_up(~isnan(valid_up));
    if ~isempty(valid_up)
        y_1(i) = prctile(valid_up, 99.95) - 30;
        y_2(i) = prctile(valid_up, 0.05) - 30;
        y_3(i) = mean(valid_up) - 30;
    else
        y_1(i) = NaN; y_2(i) = NaN; y_3(i) = NaN;
    end
    x(i) = (i-1) * step / 1950.70866;
end

% %% 观测数据
% t_check = [5.8741E+01 1.3869E+02  2.1810E+02 2.4204E+02 2.8120E+02];
% z       = [2.4882E+00 4.0190E+00  1.3130E+01 1.6057E+01 2.0226E+01] + 30;

%% 绘制结果 - 裂纹长度图
figure(1);
plot(t_check, z-30, '^', 'linewidth', 5); hold on
plot(x, y_3, 'b', 'linewidth', 5); hold on;
plot(x, y_1, 'r--', 'linewidth', 5); hold on
plot(x, y_2, 'r--', 'linewidth', 5); hold off;


%% 图表格式设置 - 裂纹长度图
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

%% 绘制结果 - 失效概率(POF)变化图
figure(2);
% 提取有效的时间步和POF值
x_pof = x(1:m-1);
pof_valid = POF_array(1:m-1);

% 绘制POF曲线
semilogy(x_pof, pof_valid, 'b-', 'linewidth', 3); hold on;

% 绘制临界POF线（10^-7）
critical_POF = 1e-7;
semilogy([x_pof(1), x_pof(end)], [critical_POF, critical_POF], 'r--', 'linewidth', 2);
hold off;

% 图表格式设置
xlabel('Flight hours/h', 'FontSize', 30);
ylabel('Probability of Failure (POF)', 'FontSize', 30);
set(get(gca, 'xlabel'), 'fontname', 'Times New Roman');
set(get(gca, 'ylabel'), 'fontname', 'Times New Roman');
set(gca, 'fontname', 'Times New Roman');
set(gca, 'FontSize', 40);

% 图例设置
legend('   POF', '   Critical POF (10^{-7})', 'FontSize', 40);
legend('boxoff')
legend('Location', 'best');

% 网格设置
grid on;
set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
set(gca, 'YScale', 'log');  % 确保使用对数坐标
