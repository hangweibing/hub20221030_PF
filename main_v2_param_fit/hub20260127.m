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

% 真实材料参数 (可选，用于对比，若没有则设为 [])
TRUE_PARAMS = [-9.7360, 2.6089, 0.1603, 0.3940]; 
true_log_theta1 = TRUE_PARAMS(1);
true_theta2 = TRUE_PARAMS(2);
true_theta3 = TRUE_PARAMS(3);
true_k2 = TRUE_PARAMS(4);

% 添加路径（可选）
% addpath surModelPackage
addpath(genpath('For_PF_260131_1'))      % 添加Param_fit及其子文件夹

% 初始化 KDE 采样
kde_sampling('preprocess');


%% ===================================================================
%% 粒子滤波算法参数配置
%% ===================================================================

% Debug模式开关
DEBUG_MODE = false;                    % 设置为 true 开启debug模式，false 关闭

% 参数拟合模式配置
fit_fix_mode = 'theta3';                    % 'k2': 固定 k2, 优化 theta; 'theta3': 固定 theta3, 优化 k2; 'none': 都不固定
k2_fixed_val = 0.394;                       % 当 fit_fix_mode='k2' 时生效
theta3_fixed_val = 0.1603;                  % 当 fit_fix_mode='theta3' 时生效

% 基本粒子滤波参数
n = 1;                                    % 状态向量的维度（每个粒子）
N = 100;                                  % 粒子数量
step = 1000;                              % 每步循环次数
v_sphere = 2;                             % 一维空间维度参数

% 正则化粒子滤波参数计算
A = (8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h = A*N^(-1/(n+4));                       % 正则化粒子滤波的平滑参数

% 预采样所有粒子的参数组合（提高效率）
SAVE_INTERVAL = 50;                       % 数据保存间隔：每隔 k 次循环保存一次数据到磁盘
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


% 内存中的粒子状态缓冲区（仅保留当前步和下一步，极大地节省 RAM）
xparticle_curr = zeros(N, 46);
% xparticle_next = zeros(N, 46);

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
t_check = [1.0000e+02 6.0591e+02 1.1118e+03 1.6177e+03 2.1236e+03 2.6295e+03 3.1355e+03 3.6414e+03 4.1473e+03 4.6532e+03];
z       = [3.2041e+01 3.2247e+01 3.2451e+01 3.2678e+01 3.2988e+01 3.3537e+01 3.5950e+01 4.0230e+01 4.4797e+01 5.0000e+01];

% t_check = [20, 60, 100, 140, 180, 220, 260, 300, 340];
% z = [13.194, 14.381, 15.281, 16.063, 16.888, 17.622, 18.341, 19.044, 19.694] + 30;

% 统计量和观测矩阵（规模较小，可保留在内存中）
% 计算大概需要的总步数，用于预分配空间
max_hours = t_check(end);
estimated_steps = ceil(max_hours * cyclesperhour / step) + 100; 

Xpf = zeros(46, estimated_steps);                   % 滤波估计值均值
% xparticle_cov = zeros(46, 46, 1000);     % 协方差矩阵
a_upper = zeros(1, estimated_steps);
a_lower = zeros(1, estimated_steps);
POF_array = zeros(estimated_steps, 1);              % 失效概率数组
% D = zeros(46, 1000);                     % 协方差对角线（RPF残余）


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

    % --- 物理一致性修正：强制执行参数固定逻辑 ---
    if strcmp(fit_fix_mode, 'k2')
        k2 = k2_fixed_val;       % 强制所有粒子使用固定的 k2
    elseif strcmp(fit_fix_mode, 'theta3')
        theta3 = theta3_fixed_val; % 强制所有粒子使用固定的 theta3
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

% 打印状态以验证
if strcmp(fit_fix_mode, 'k2')
    fprintf('  [物理约束] 已强制固定所有粒子的 k2 = %.4f\n', k2_fixed_val);
elseif strcmp(fit_fix_mode, 'theta3')
    fprintf('  [物理约束] 已强制固定所有粒子的 theta3 = %.4f\n', theta3_fixed_val);
end

% 初始化粒子权重
current_weight = 1/N * ones(N, 1);

% 将初始时刻数据保存到缓冲区（第1步）
buffer_step_counter = 1;
% buffer_xparticle(:, :, 1) = xparticle_curr;
buffer_upcrackparticles(:, 1) = xparticle_curr(:, 42);

% --- [修复] 初始化初始时刻的统计量 ---
valid_mask_init = ~isnan(xparticle_curr(:, 1));
if any(valid_mask_init)
    Xpf(:, 1) = mean(xparticle_curr(valid_mask_init, :))';
    a_upper(1) = prctile(xparticle_curr(valid_mask_init, 42), 99.9);
    a_lower(1) = prctile(xparticle_curr(valid_mask_init, 42), 0.1);
end
% -------------------------------------

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
%% 参数拟合初始化（用于观测步）
%% ===================================================================
% 加载参数拟合所需的 GP 模型和人口数据
load('For_PF_260131_1/AM-TC4-GRO_260123_merged.mat');
load('For_PF_260131_1/parameter_gp_AM-TC4-GRO_combine.mat');
p_model = parameter_gp;  % 创建副本，与 demo 代码保持一致
p_model.spectra = spectra; % 使用当前主程序的载荷谱

% 初始化优化参数
options_fmincon_1 = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'SpecifyObjectiveGradient',false, ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 2000, ...
    'Display', 'off');

options_fmincon_2 = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'SpecifyObjectiveGradient',false, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 1000, ...
    'Display', 'off');

p_model.fmincon_option1 = options_fmincon_1;
p_model.fmincon_option2 = options_fmincon_2;

% 设置初始参数
p_model.PARA0 = pop_now{1,6};
p_model.PARA0 = (p_model.PARA0(:))'; % 强制转换为行向量

% 参数固定逻辑设置
if strcmp(fit_fix_mode, 'k2')
    p_model.fix_k2 = true; 
    p_model.k2_fixed_val = k2_fixed_val;
    p_model.fix_theta3 = false;
    fprintf('[参数固定] 已设置固定 k2 = %.4f\n', k2_fixed_val);
elseif strcmp(fit_fix_mode, 'theta3')
    p_model.fix_k2 = false; 
    p_model.fix_theta3 = true;
    p_model.theta3_fixed_val = theta3_fixed_val;
    fprintf('[参数固定] 已设置固定 theta3 = %.4f\n', theta3_fixed_val);
else
    p_model.fix_k2 = false;
    p_model.fix_theta3 = false;
    fprintf('[参数固定] 未固定任何参数，所有参数将参与优化\n');
end

% 初始化边界
theta_init = p_model.PARA0(3:end);
theta_gene = theta_init(2:end-1);
f_L = zeros(p_model.numGenes, 1);
f_U = zeros(p_model.numGenes, 1);
for pp = 1:p_model.numGenes
    if theta_gene(pp) > 0
        f_L(pp) = 0.1; f_U(pp) = 10;
    elseif theta_gene(pp) < 0
        f_L(pp) = -10; f_U(pp) = -0.1;
    else
        f_L(pp) = 0; f_U(pp) = 0;
    end
end
p_model.LB_orig = [-15; f_L; 1; 0; 0];
p_model.UB_orig = [-5; f_U; 10; 500; 500];

% 处理 evalstr
evalstr = p_model.evalstr1;
if iscell(evalstr), evalstr = evalstr{1}; end % 确保是字符串而非 cell
evalstr = regexprep(evalstr, 'c(\d+)', 'Const_pair_now($1)');
evalstr = regexprep(evalstr, 'x(\d+)', 'xtrain(:,$1)');
p_model.evalstr2 = evalstr;

% 其他必要的参数设置
p_model.parameter_K = p_model.PARA0(1:2);
p_model.theta = theta_init;
p_model.theta_end_limit = [-100, 100]; % 补齐缺失的约束范围变量
p_model.ytrain = 0.7; 
p_model.cyclesperhour = cyclesperhour;

% 封装仿真相关参数，以便传递给 pred_a_N
p_model.spectrum_factor = 60;
p_model.ref_load = 100;
p_model.W = 60;
p_model.B = 5;
p_model.step_size = step;

% 将 p_model 赋值回 parameter_gp 以便后续使用
parameter_gp = p_model;

% (平均应力计算已移除，改为直接在 a2aFunc 中循环处理载荷谱)

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
prev_PARA0_fit = []; % 用于存储拟合结果历史
while (m-1)*step/cyclesperhour <= t_check(end)
    iter_tic = tic;  % 初始化当前步耗时计时器
    xparticle_prev = xparticle_curr;
    %% 时间步数据准备 (载荷谱段提取)
    % 提取当前 step (1000个循环) 的载荷段
    idx_start = 2 * (m-2) * step + 1;
    idx_end = 2 * (m-1) * step;
    loads_segment = spectra(idx_start : idx_end);

    %% ===================================================================
    %% 粒子预测步骤（对每个粒子进行状态更新）
    %% ===================================================================

    % 初始化当前时间步临时存储
    particles_deltaK_max_temp = zeros(N, 1);
    x_next_temp = zeros(N, 46);
    splitted_update = zeros(1, N);

    for (i = 1:N)
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

            %% 调用预测模型更新粒子状态 (传入载荷谱段以支持循环计算)
            [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, deltaKSet] = ...
                a2aNew(yRegSet, zRegSet, loads_segment, m_name, ...
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

    if any(valid_mask)
        Xpf(:, m) = (mean(xparticle_curr(valid_mask, :)))';
        a_upper(m) = prctile(xparticle_curr(valid_mask, 42), 99.9);
        a_lower(m) = prctile(xparticle_curr(valid_mask, 42), 0.1);
        % xparticle_cov(:, :, m) = cov(xparticle_curr(valid_mask, :));
    else
        Xpf(:, m) = Xpf(:, m-1); % 如果全部失效，保持不变
        a_upper(m) = a_upper(m-1);
        a_lower(m) = a_lower(m-1);
    end

    %% 计算POF（！！！！！注意K的单位换算！！！！！）
    try
        valid_K = particles_deltaK_max_temp(~isnan(particles_deltaK_max_temp));
        if ~isempty(valid_K)
            % 这里假设 K 已通过 a2aFunc 内部处理为有效的 Kmax (或 ΔK/1-R)
            POF_array(m) = calculatePOF(valid_K / sqrt(1000), mu_Kc, std_Kc);
        else
            POF_array(m) = POF_array(m-1);
        end
    catch ME
        POF_array(m) = 0;
    end

    %% 观测更新步骤（当到达观测时刻时）
    if ((m-1)*step/cyclesperhour < t_check(j)) && (m*step/cyclesperhour >= t_check(j))
        fprintf('  [观测更新] 时间: %.2f h, 观测值: %.2f mm\n', t_check(j), z(j));

        % --- [Multi-Weight] 1. 参数拟合 (得到参数的"观测值") ---
        % 初始化多维度权重矩阵 [N x 5] (z_up, log_theta1, theta2, theta3, k2)
        multi_weights = ones(N, 5) / N; 
        
        % 在第一次观测时，保存权重最大粒子的几何信息作为拟合起始点
        if j == 1
            % 计算裂纹长度权重（仅用于找最佳粒子）
            zPred_curr = xparticle_curr(:, 42);
            weights_temp = zeros(N, 1);
            for i = 1:N
                if isnan(zPred_curr(i))
                    weights_temp(i) = 1e-99;
                else
                    res = z(j) - zPred_curr(i);
                    weights_temp(i) = (1/sqrt(2*pi*R)) * exp(-0.5*(res^2)/R) + 1e-99;
                end
            end
            % 找到权重最大的粒子
            [~, best_idx] = max(weights_temp);
            % 保存该粒子的几何信息 (y坐标集 + z坐标集)
            parameter_gp.init_geometry = struct();
            parameter_gp.init_geometry.yRegSet = xparticle_curr(best_idx, 1:21);
            parameter_gp.init_geometry.zRegSet = xparticle_curr(best_idx, 22:42);
            fprintf('  [初始几何] 已保存粒子 %d 的几何信息作为拟合起始点 (z_end=%.4f)\n', ...
                best_idx, xparticle_curr(best_idx, 42));
        end
        
        fit_triggered = false;
        if j > 1
            fprintf('  [参数拟合] 正在根据前 %d 次观测数据进行多起始点拟合...\n', j);
            
            % 准备拟合所需的历史观测数据 (1 到 j)
            parameter_gp.data_a_N = [t_check(1:j)', z(1:j)'];
            parameter_gp.num_of_data_a_N = size(parameter_gp.data_a_N, 1);
            
            % 执行多起始点优化
            N_starts_opt = 10;
            [~, delta_kth_fit, kc_fit, theta_fit, ~, ~, pass_idx] = ...
                loss_cal_optimize_for_PF_multistart(parameter_gp, pop_now, N_starts_opt);
            
            if pass_idx == 1
                % 注意: theta_fit = [log_theta1, theta3, theta2]
                fprintf('  [参数拟合] 拟合成功！结果：f1=%.4f, f2(theta3)=%.4f, f3(theta2)=%.4f, kc=%.4f\n', ...
                    theta_fit(1), theta_fit(2), theta_fit(3), kc_fit);
                
                % 保存拟合结果供下一次拟合参考
                new_fit = [delta_kth_fit, kc_fit, theta_fit'];
                prev_PARA0_fit = [prev_PARA0_fit; new_fit];
                
                % 设置参数观测值
                % theta_fit 顺序 (PARA0): [log_theta1, theta3, theta2]
                % 粒子状态顺序: [log_theta1, theta2, theta3, k2]
                % 需要重新映射以匹配粒子状态顺序
                param_obs = [theta_fit(1), theta_fit(3), theta_fit(2), kc_fit]; % [log_theta1, theta2, theta3, kc]
                fit_triggered = true;
            else
                fprintf('  [参数拟合] 拟合未通过约束检查。\n');
            end
        end

        % --- [Multi-Weight] 2. 计算各维度似然权重 ---
        % 2.1 裂纹长度权重 (z)
        zPred_curr = xparticle_curr(:, 42);
        for i = 1:N
            if isnan(zPred_curr(i))
                multi_weights(i, 1) = 1e-99;
            else
                res_z = z(j) - zPred_curr(i);
                multi_weights(i, 1) = (1/sqrt(2*pi*R)) * exp(-0.5*(res_z^2)/R) + 1e-99;
            end
        end

        % 2.2 模型参数权重 (仅在拟合成功时计算)
        if fit_triggered
            for d = 2:5
                obs_val = param_obs(d-1);
                pred_vals = xparticle_curr(:, 43 + (d-2)); % 43: log_theta1, 44: theta2, 45: theta3, 46: k2
                for i = 1:N
                    if isnan(pred_vals(i))
                        multi_weights(i, d) = 1e-99;
                    else
                        res_p = obs_val - pred_vals(i);
                        % 使用与裂纹长度相同的噪声或独立定义
                        multi_weights(i, d) = (1/sqrt(2*pi*R)) * exp(-0.5*(res_p^2)/R) + 1e-99;
                    end
                end
            end
        end

        % --- [Multi-Weight] 3. 可视化似然权重分布 ---
        if ~exist('param_obs', 'var'), param_obs = []; end
        xpart_for_plot = [xparticle_curr(:, 42), xparticle_curr(:, 43:46)];
        plot_multi_weight_distribution(xpart_for_plot, multi_weights, j, t_check(j), z(j), param_obs, TRUE_PARAMS);

        % --- [Multi-Weight] 4. 联合权重重采样 ---
        resample_opts = struct();
        resample_opts.weight_mode = 'product';
        resample_opts.kernel_scale = 0.5;
        
        % 注意：重采样需在绘图之后进行
        xparticle_curr = resample_particles_joint(xparticle_curr, multi_weights, resample_opts);
        
        current_weight = 1/N * ones(N, 1);

        % 重采样后重新计算统计量并更新
        valid_mask_res = ~isnan(xparticle_curr(:, 1));
        if any(valid_mask_res)
            Xpf(:, m) = mean(xparticle_curr(valid_mask_res, :))';
            a_upper(m) = prctile(xparticle_curr(valid_mask_res, 42), 99.9);
            a_lower(m) = prctile(xparticle_curr(valid_mask_res, 42), 0.1);
        end

        % --- [新增] 实时绘制裂纹长度图并保存 ---
        fig_rt = figure('Name', sprintf('Crack Length Realtime - Obs %d', j), 'Visible', 'off', 'Position', [100, 100, 1000, 600], 'Color', 'w');
        cur_time_axis = (0:m-1) * step / cyclesperhour;
        plot(t_check(1:j), z(1:j), '^', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k'); hold on;
        plot(cur_time_axis, Xpf(42, 1:m), 'b', 'LineWidth', 3); hold on;
        plot(cur_time_axis, a_upper(1:m), 'r--', 'LineWidth', 2); hold on;
        plot(cur_time_axis, a_lower(1:m), 'r--', 'LineWidth', 2); hold off;
        xlabel('Flight hours/h', 'FontSize', 12); ylabel('Crack length/mm', 'FontSize', 12);
        legend('Experimental value', 'Prediction mean', '99.9% bounds', 'Location', 'best');
        title(sprintf('Crack Growth Real-time Prediction (Obs %d, Time %.2f h)', j, t_check(j)));
        grid on; set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
        
        save_path_rt = 'MultiWeight_Dist_Results';
        if ~exist(save_path_rt, 'dir'), mkdir(save_path_rt); end
        saveas(fig_rt, fullfile(save_path_rt, sprintf('CrackLength_Obs_%d.png', j)));
        close(fig_rt);

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
    if buffer_step_counter >= SAVE_INTERVAL || (m*step/cyclesperhour > t_check(end))
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

    disp(['已完成' num2str((m-1)*step/cyclesperhour) '小时，进行了' num2str(j-1) '次观测', ...
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
    x(i) = (i-1) * step / cyclesperhour;
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
