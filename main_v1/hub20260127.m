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
SAVE_INTERVAL = 10;                       % 数据保存间隔：每隔 k 次循环保存一次数据到磁盘
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
% 粒子滤波数据存储配置
cyclesperhour = 1950.70866;              % 每小时循环次数
mu_Kc = 33.4 ;                           % 断裂韧性均值 (MPa√m)
std_Kc = 3.34 ;                          % 断裂韧性标准差 (MPa√m)

% 使用 matfile 进行批量存储（先在内存中累积，再批量写入磁盘）
% 避免频繁 I/O 操作，提升性能
history_file = 'pf_simulation_results.mat';
if exist(history_file, 'file'), delete(history_file); end
mfile = matfile(history_file, 'Writable', true);

% 在磁盘上预分配空间
% mfile.xparticle = zeros(N, 46, 2000);
mfile.upcrackparticles = zeros(N, 2000);
% mfile.log_theta1_particles = zeros(N, 2000);
% mfile.theta2_particles = zeros(N, 2000);
% mfile.theta3_particles = zeros(N, 2000);
% mfile.k2_particles = zeros(N, 2000);
% mfile.weight = zeros(N, 2000);

% 内存缓冲区：用于累积 SAVE_INTERVAL 步的数据，减少 I/O 频率
% buffer_xparticle = zeros(N, 46, SAVE_INTERVAL);
buffer_upcrackparticles = zeros(N, SAVE_INTERVAL);
% buffer_log_theta1_particles = zeros(N, SAVE_INTERVAL);
% buffer_theta2_particles = zeros(N, SAVE_INTERVAL);
% buffer_theta3_particles = zeros(N, SAVE_INTERVAL);
% buffer_k2_particles = zeros(N, SAVE_INTERVAL);
% buffer_weight = zeros(N, SAVE_INTERVAL);
buffer_step_counter = 0;  % 缓冲区当前累积的步数

% 统计量和观测矩阵（规模较小，可保留在内存中）
% Xpf = zeros(46, 1000);                   % 滤波估计值均值
% xparticle_cov = zeros(46, 46, 1000);     % 协方差矩阵
POF_array = zeros(1000, 1);              % 失效概率数组
% D = zeros(46, 1000);                     % 协方差对角线（RPF残余）

% 内存中的粒子状态缓冲区（仅保留当前步和下一步，极大地节省 RAM）
xparticle_curr = zeros(N, 46);
% xparticle_next = zeros(N, 46);

%% ===================================================================
%% 观测数据配置
%% ===================================================================

% 铝合金（真实数据）
% t_check=[5.8741E+01 9.9534E+01 1.3869E+02 171.873 1.8275E+02 203.419 2.1810E+02 2.4204E+02 2.8120E+02 3.1546E+02 3.4375E+02 3.8182E+02];
% z      =[2.4882E+00 2.9894E+00 4.0190E+00 5.46751 7.2212E+00 10.4231 1.3130E+01 1.6057E+01 2.0226E+01 2.1715E+01 2.2903E+01 2.3452E+01]+30;

% demo测试数据（钛合金）
t_check = [1.0000e+02 6.0591e+02 1.1118e+03 1.6177e+03 2.1236e+03 2.6295e+03 3.1355e+03 3.6414e+03 4.1473e+03 4.6532e+03];
z       = [3.2041e+01 3.2247e+01 3.2451e+01 3.2678e+01 3.2988e+01 3.3537e+01 3.5950e+01 4.0230e+01 4.4797e+01 5.0000e+01];

% 试验数据（铝合金）
% t_check = [20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340];
% z = [14.61, 14.85, 15.24, 15.43, 15.95, 16.22, 16.57, 16.88, 17.17, 17.69, 17.83, 18.42, 18.77, 19.01, 19.35, 19.8, 20.34];
% z = [13.194, 13.656, 14.381, 14.869, 15.281, 15.709, 16.063, 16.553, 16.888, 17.303, 17.622, 17.869, 18.341, 18.791, 19.044, 19.394, 19.694, 20.038];
t_check = [20, 60, 100, 140, 180, 220, 260, 300, 340];
z_half = [13.194, 14.381, 15.281, 16.063, 16.888, 17.622, 18.341, 19.044, 19.694];
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

    %% 粒子状态存储 (存入当前步缓冲区)
    % 完整状态向量：[y坐标集(21), z坐标集(21), log_theta1_, theta2, theta3, k2]
    xparticle_curr(i, :) = [yIniRegSet, zIniRegSet, log_theta1_, theta2, theta3, k2];
end

% 初始化粒子权重
current_weight = 1/N * ones(N, 1);

% 将初始时刻数据保存到缓冲区（第1步）
buffer_step_counter = 1;
% buffer_xparticle(:, :, 1) = xparticle_curr;
buffer_upcrackparticles(:, 1) = xparticle_curr(:, 42);
% buffer_log_theta1_particles(:, 1) = xparticle_curr(:, 43);
% buffer_theta2_particles(:, 1) = xparticle_curr(:, 44);
% buffer_theta3_particles(:, 1) = xparticle_curr(:, 45);
% buffer_k2_particles(:, 1) = xparticle_curr(:, 46);
% buffer_weight(:, 1) = current_weight;

%% ===================================================================
%% 初始时刻统计量计算
%% ===================================================================


% 初始化粒子坐标信息存储（cell数组，1×N，每列保存一个粒子的坐标历史）
particles_coordinates = cell(1, N);

% % 计算初始滤波估计（各参数的均值）
% Xpf(:, 1) = (mean(xparticle_curr))';
%
% % 计算初始协方差矩阵
% xparticle_cov(:, :, 1) = cov(xparticle_curr);

%% ===================================================================
%% 数据加载和预处理
%% ===================================================================
% 加载ASTM标准载荷谱数据
load('AsteixSpectraData.mat');
load('pod_models.mat'); % 加载POD模型数据

% 设置载荷谱重复次数，确保仿真时长覆盖观测点
spectra_repeat_count = 80;
spectra = repmat(reshape(spectra, 1, []), 1, spectra_repeat_count);

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
    % 获取上一时刻的所有粒子状态 (此时 xparticle_curr 实际上是 m-1 步的结果)
    xparticle_prev = xparticle_curr;

    % 获取当前时间段的平均应力增量
    aver_delta_sigma = aver_delta_sigma_set(m-1);

    %% ===================================================================
    %% 粒子预测步骤（对每个粒子进行状态更新）
    %% ===================================================================

    % 初始化当前时间步临时存储
    particles_deltaK_max_temp = zeros(N, 1);
    x_next_temp = zeros(N, 46);
    splitted_update = zeros(1, N);

    parfor (i = 1:N)
        try
            %% 粒子级变量初始化
            curUinput = {};
            curAverInput = {};
            SPLITTED = SPLITTE_temp(i);  % 获取上一时刻的分裂状态

            %% 为每个并行线程设置独立的随机数流
            stream = RandStream('mt19937ar', 'Seed', SIM_SEED + i + m*N);

            %% 粒子状态提取
            xparticlei = xparticle_prev(i, :);

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
                    xparticlei = xparticle_prev(tmp_sel, :);
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
            % 记录坐标历史 (仅在 Debug 模式开启时记录，避免非必要的内存消耗)
            if DEBUG_MODE
                particles_coordinates{i} = [particles_coordinates{i}; [yRegSet, zRegSet]];
            end

            log_theta1_ = xparticlei(43);
            theta2 = xparticlei(44);
            theta3 = xparticlei(45);
            k2 = xparticlei(46);

            %% 调用预测模型更新粒子状态
            [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, deltaKSet] = ...
                a2aNew(yRegSet, zRegSet, aver_delta_sigma, aver_R_set(m-1), aver_Smax_set(m-1), m_name, ...
                curUinput, curAverInput, log_theta1_, theta2, theta3, k2, step, testErrSet, i);

            particles_deltaK_max_temp(i) = deltaKSet(end);
            x_next_temp(i, :) = real([yRegSet, zRegSet, log_theta1_, theta2, theta3, k2]);
            splitted_update(i) = SPLITTED;

        catch ME
            x_next_temp(i, :) = NaN;
            particles_deltaK_max_temp(i) = NaN;
            splitted_update(i) = 0;
        end
    end

    % 更新缓冲区
    xparticle_curr = x_next_temp;
    SPLITTE_temp = splitted_update;

    % 计算粒子滤波统计量 (排除 NaN 无效粒子)
    valid_mask = ~isnan(xparticle_curr(:, 1));
    % 默认本步权重（观测前）为均匀分布
    current_weight = 1/N * ones(N, 1);

    %     if any(valid_mask)
    %         Xpf(:, m) = (mean(xparticle_curr(valid_mask, :)))';
    %         xparticle_cov(:, :, m) = cov(xparticle_curr(valid_mask, :));
    %     else
    %         Xpf(:, m) = Xpf(:, m-1); % 如果全部失效，保持不变
    %     end

    %% 计算POF（！！！！！注意K的单位换算！！！！！）
    try
        valid_K = particles_deltaK_max_temp(~isnan(particles_deltaK_max_temp));
        if ~isempty(valid_K)
            particles_K_max_valid = valid_K / (1 - aver_R_set(m-1));
            POF_array(m) = calculatePOF(particles_K_max_valid / sqrt(1000), mu_Kc, std_Kc);
        else
            POF_array(m) = POF_array(m-1);
        end
    catch ME
        POF_array(m) = 0;
    end

    %% 观测更新步骤（当到达观测时刻时）
    if ((m-1)*step/1950.70866 < t_check(j)) && (m*step/1950.70866 >= t_check(j))
        %% 计算似然权重
        zPred_curr = xparticle_curr(:, 42);
        weights_temp = zeros(N, 1);
        for i = 1:N
            if isnan(zPred_curr(i))
                weights_temp(i) = 1e-99;
                continue;
            end
            res = z(j) - zPred_curr(i);
            weights_temp(i) = (1/sqrt(2*pi*R)) * exp(-0.5*(res^2)/R) + 1e-99;
        end

        %% 归一化权重
        current_weight = weights_temp ./ sum(weights_temp);

        %         %% 更新状态估计
        %         Xpf(:, m) = (xparticle_curr' * current_weight);
        %
        %         %% 更新协方差
        %         diff = xparticle_curr - Xpf(:, m)';
        %         xparticle_cov(:, :, m) = (diff' .* current_weight') * diff;

        % %% 正则化重采样
        % for i = 1:44
        %     D(i, m) = sqrt(xparticle_cov(i, i, m));
        % end
        outindex = randomr(current_weight);
        xparticle_curr = xparticle_curr(outindex, :);
        current_weight = 1/N * ones(N, 1); % 重采样后权重重置

        % 重采样后清空粒子坐标历史
        particles_coordinates = cell(1, N);
        j = j + 1;
    end

    %% ===================================================================
    %% 数据累积与批量保存（减少 I/O 消耗）
    %% ===================================================================

    % 每次循环都将数据累积到内存缓冲区
    buffer_step_counter = buffer_step_counter + 1;
    %     buffer_xparticle(:, :, buffer_step_counter) = xparticle_curr;
    buffer_upcrackparticles(:, buffer_step_counter) = xparticle_curr(:, 42);
    %     buffer_log_theta1_particles(:, buffer_step_counter) = xparticle_curr(:, 43);
    %     buffer_theta2_particles(:, buffer_step_counter) = xparticle_curr(:, 44);
    %     buffer_theta3_particles(:, buffer_step_counter) = xparticle_curr(:, 45);
    %     buffer_k2_particles(:, buffer_step_counter) = xparticle_curr(:, 46);
    %     buffer_weight(:, buffer_step_counter) = current_weight;

    % 每隔 SAVE_INTERVAL 步批量保存一次数据，或者在最后一步强制保存
    if buffer_step_counter >= SAVE_INTERVAL || (m*step/1950.70866 > t_check(end))
        % 计算本次要保存的时间步范围
        start_step = m - buffer_step_counter + 1;
        end_step = m;

        % 批量保存粒子状态到磁盘 (matfile)
        %         mfile.xparticle(:, :, start_step:end_step) = buffer_xparticle(:, :, 1:buffer_step_counter);
        mfile.upcrackparticles(:, start_step:end_step) = buffer_upcrackparticles(:, 1:buffer_step_counter);
        %         mfile.log_theta1_particles(:, start_step:end_step) = buffer_log_theta1_particles(:, 1:buffer_step_counter);
        %         mfile.theta2_particles(:, start_step:end_step) = buffer_theta2_particles(:, 1:buffer_step_counter);
        %         mfile.theta3_particles(:, start_step:end_step) = buffer_theta3_particles(:, 1:buffer_step_counter);
        %         mfile.k2_particles(:, start_step:end_step) = buffer_k2_particles(:, 1:buffer_step_counter);
        %         mfile.weight(:, start_step:end_step) = buffer_weight(:, 1:buffer_step_counter);

        % 如果开启了 Debug 模式，保存坐标历史到 mat 文件
        if DEBUG_MODE
            save('debug_particles_coordinates.mat', 'particles_coordinates', 'm');
        end

        fprintf('  [数据保存] 已批量保存时间步 %d-%d 的粒子数据到磁盘 (%d 步)\n', ...
            start_step, end_step, buffer_step_counter);

        % 重置缓冲区计数器
        buffer_step_counter = 0;
    end

    m = m + 1;  % 时间步递增

    % 计算耗时
    iter_time = toc(iter_tic);
    total_time = toc(total_tic);

    % 计算当前步的平均裂纹长度（仅用于显示进度）
    mean_z = mean(xparticle_curr(:, 42), 'omitnan');

    disp(['已完成' num2str((m-1)*step/1950.70866) '小时，进行了' num2str(j-1) '次观测', ...
        '，当前z均值：' num2str(mean_z, '%.4f'), ...
        '，当前步耗时：' num2str(iter_time, '%.2f') 's', ...
        '，累计总耗时：' num2str(total_time, '%.2f') 's']);
end

%% ===================================================================
%% 结果后处理和可视化
%% ===================================================================
close all
%% 数据处理
clear x y_0 y_1 y_2 y_3 y_33 y21 PoF3

% 从磁盘提取历史数据用于后处理
upcrack_history = mfile.upcrackparticles;

% 计算统计量
for i = 1:m-1
    valid_up = upcrack_history(:, i);
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

% 保存图像
saveas(gcf, 'result_crack_length.png');
saveas(gcf, 'result_crack_length.fig');

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

% 保存图像
saveas(gcf, 'result_POF.png');
saveas(gcf, 'result_POF.fig');
