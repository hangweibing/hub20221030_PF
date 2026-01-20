%% ===================================================================
%% 程序初始化和环境设置
%% ===================================================================

% 清理工作空间
clear; close all; clc; tic; format long;

%% ===================================================================
%% 随机数种子设置（确保结果可重复）
%% ===================================================================
SIM_SEED = 2023;  % 您可以修改此数值来获得不同的随机序列
rng(SIM_SEED);    % 设置全局随机数种子

% 添加路径（可选）
% addpath surModelPackage

%% ===================================================================
%% 粒子滤波算法参数配置
%% ===================================================================

% 基本粒子滤波参数
n = 1;                                    % 状态向量的维度（每个粒子）
N = 20000;                                  % 粒子数量
v_sphere = 2;                             % 一维空间维度参数

% 正则化粒子滤波参数计算
A = (8/v_sphere*(n+4)*(2*sqrt(pi))^n)^(1/(n+4));
h = A*N^(-1/(n+4));                       % 正则化粒子滤波的平滑参数

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
logDparticles = zeros(N, 1000);          % logD参数粒子（NASGRO模型）
Aparticles = zeros(N, 1000);             % A参数粒子（NASGRO模型）
delta_kthrparticles = zeros(N, 1000);     % delta_kthr参数粒子（NASGRO模型）
pparticles = zeros(N, 1000);             % p参数粒子（NASGRO模型）
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

% 实际使用的观测数据（简化版）
t_check = [5.8741E+01 9.9534E+01 1.3869E+02  1.8275E+02 2.4204E+02 2.8120E+02];
z       = [2.4882E+00 2.9894E+00 4.0190E+00  7.2212E+00 1.6057E+01 2.0226E+01] + 30;

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
    %% 参数采样（NASGRO(H-S)模型四参数）
    % logD: 均匀分布 
    logD = unifrnd(-10.5, -9.00, 1, 1);
    
    % A: 均匀分布 
    A = unifrnd(110, 160, 1, 1);
    
    % delta_kthr: 均匀分布
    delta_kthr = unifrnd(0.1, 1.5, 1, 1);
    
    % p: 正态分布 
    % p = normrnd(2.12, 0.04, 1, 1);
    p = unifrnd(1.80, 2.18, 1, 1);
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
    % 完整状态向量：[y坐标集(21), z坐标集(21), logD, A, delta_kthr, p]
    xparticle(i, :, 1) = [yIniRegSet, zIniRegSet, logD, A, delta_kthr, p];
end

%% ===================================================================
%% 初始时刻统计量计算
%% ===================================================================

% 提取各参数的粒子分布
upcrackparticles(:, 1) = xparticle(:, 42, 1);      % 上表面裂纹长度
logDparticles(:, 1) = xparticle(:, 43, 1);         % logD参数（NASGRO模型）
Aparticles(:, 1) = xparticle(:, 44, 1);            % A参数（NASGRO模型）
delta_kthrparticles(:, 1) = xparticle(:, 45, 1);  % delta_kthr参数（NASGRO模型）
pparticles(:, 1) = xparticle(:, 46, 1);            % p参数（NASGRO模型）

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

% Debug模式开关
DEBUG_MODE = true;                    % 设置为 true 开启debug模式，false 关闭

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
        %% 粒子级变量初始化
        curUinput = {};
        curAverInput = {};
        SPLITTED = SPLITTE_temp(i);  % 获取上一时刻的分裂状态
        
        %% 为每个并行线程设置独立的随机数流（parfor线程安全）
        % 原方案: 直接使用全局randi()，导致线程间竞态条件和不可重现结果
        % 问题: parfor中所有线程共享全局随机数生成器状态
        % 解决: 为每个线程创建独立随机流，避免线程间干扰
        % 使用全局种子、粒子索引和时间步作为复合种子，确保完全可重现性
        stream = RandStream('mt19937ar', 'Seed', SIM_SEED + i + m*N);

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
        %% 模型选择
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

        % 保存当前粒子的坐标信息（第m次更新前的[y,z]坐标）
        particles_coordinates{i} = [particles_coordinates{i}; [yRegSet, zRegSet]];

        % NASGRO(H-S)模型参数
        logD = xparticlei(43);           % NASGRO模型参数logD
        A = xparticlei(44);              % NASGRO模型参数A
        delta_kthr = xparticlei(45);     % NASGRO模型参数delta_kthr
        p = xparticlei(46);              % NASGRO模型参数p

        % 可选：绘制裂纹几何形状
        % figure; plot_geometry_20; plot(zRegSet, yRegSet); axis equal;

        %% 调用预测模型更新粒子状态
        [yRegSet, zRegSet, SPLITTED, logD, A, delta_kthr, p, deltaKSet] = ...
            a2aNew(yRegSet, zRegSet, aver_delta_sigma, aver_R_set(m-1), m_name, ...
                   curUinput, curAverInput, logD, A, delta_kthr, p, step, testErrSet, i);

        %% 存储当前粒子的100%分位数应力强度因子
        particles_deltaK_max_temp(i) = max(deltaKSet);

        %% 更新粒子状态（parfor兼容：直接确保实数）
        % 原方案: xparticle(i, :, m) = [...]; xparticle = real(xparticle);
        % 问题: 第二行违反parfor切片规则，读取整个数组导致竞态条件
        % 解决: 在赋值时直接取实部，避免全局数组操作
        xparticle(i, :, m) = real([yRegSet, zRegSet, logD, A, delta_kthr, p]);
        SPLITTE_temp(i) = SPLITTED;      % 更新分裂状态
    end

    % Debug模式：统一保存所有粒子的坐标信息
    if DEBUG_MODE
        save('debug_particles_coordinates.mat', 'particles_coordinates', 'm');
    end

    % 计算粒子滤波统计量
    weight(:, m) = 1/N * ones(N, 1);             % 均匀权重初始化
    Xpf(:, m) = (mean(xparticle(:, :, m)))';     % 状态均值估计
    xparticle_cov(:, :, m) = cov(xparticle(:, :, m)); % 状态协方差
    
    %% ===================================================================
    %% 计算POF（！！！！！注意K的单位换算！！！！！）
    %% ===================================================================
    try
        % deltaK_max 计算出 K_max
        particles_K_max = particles_deltaK_max_temp / (1 - aver_R_set(m-1));
        POF_array(m) = calculatePOF(particles_K_max / sqrt(1000), mu_Kc, std_Kc);
    catch ME
        warning('POF计算失败 (时间步 %d): %s', m, ME.message);
        POF_array(m) = 0;  % 失败时设为0
    end


    %% 提取关键参数的历史记录
    upcrackparticles(:, m) = xparticle(:, 42, m);      % 上表面裂纹长度
    logDparticles(:, m) = xparticle(:, 43, m);         % logD参数（NASGRO模型）
    Aparticles(:, m) = xparticle(:, 44, m);            % A参数（NASGRO模型）
    delta_kthrparticles(:, m) = xparticle(:, 45, m);  % delta_kthr参数（NASGRO模型）
    pparticles(:, m) = xparticle(:, 46, m);           % p参数（NASGRO模型）
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
    y_1(i) = prctile(upcrackparticles(:, i), 99.95) - 30;  % 99.95%分位数
    y_2(i) = prctile(upcrackparticles(:, i), 0.05) - 30;   % 0.05%分位数
    y_3(i) = Xpf(42, i) - 30;                              % 均值估计
    x(i) = (i-1) * step / 1950.70866;                      % 时间轴
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
