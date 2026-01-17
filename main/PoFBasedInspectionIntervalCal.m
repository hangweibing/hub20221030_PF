%% ===================================================================
%% 函数名称：PoFBasedInspectionIntervalCal
%% 功能描述：基于失效概率(POF)的检查间隔计算
%% ===================================================================
%% 理论依据：
%%   根据航空安全标准FAR-25.303，失效概率POF应小于临界值10^-7
%%   通过粒子滤波模拟裂纹扩展，计算达到临界POF所需的时间间隔
%%
%% 输入参数：
%%   xparticle: 粒子状态矩阵（从主程序传入）
%%   N: 粒子数量
%%   n_nodes: 裂纹前沿节点数量
%%   centers: 裂纹中心坐标
%%   thetas: 角度参数
%%   aver_delta_sigma_set: 应力范围序列
%%   Uinput_integrated: 完整裂纹代理模型输入
%%   Uinput_splitted_1/2: 分离裂纹代理模型输入
%%   averInput_integrated: 完整裂纹代理模型平均输入
%%   averInput_splitted_1/2: 分离裂纹代理模型平均输入
%%   step: 时间步长
%%   testErrSet: 测试误差设置
%%
%% 输出参数：
%%   InspectionInterval: 检查间隔时间（小时）
%%   PoFSet: 失效概率序列
%%
%% 作者：自动生成
%% 日期：2026-01-17
%% ===================================================================

%% ===================================================================
%% 步骤1: 初始化粒子状态用于检查间隔计算
%% ===================================================================
xparticleForCalII=xparticle;
% 生成初始裂纹尺寸分布（魏布尔分布 + 偏移）
a_ini_set_forIICal = (wblrnd(0.0057, 1.1739, 1, N) + 0.0446369) .* 25.4;

% 为每个粒子生成初始裂纹形状
for i = 1:N
    a = a_ini_set_forIICal(i);
    c = a;

    % 生成椭圆形裂纹前沿坐标
    for j = 1:n_nodes
        yIniRegSet(j) = centers(1) + c * sin(thetas(j));
        zIniRegSet(j) = centers(2) + a * cos(thetas(j));
    end

    % 更新粒子状态（y坐标在前，z坐标在后）
    xparticleForCalII(i, 1:2*n_nodes) = [yIniRegSet, zIniRegSet];
end

%% ===================================================================
%% 步骤2: 计算初始失效概率（t=0时刻）
%% ===================================================================

% 初始化分离标志
SPLITTE_temp = zeros(1, N);

% 复制粒子状态
xparticlek_1 = xparticleForCalII;

% 获取第一个应力范围
aver_delta_sigma = aver_delta_sigma_set(1);

% 并行计算每个粒子的应力强度因子增量
for i = 1:N
    % 初始化代理模型输入
    curUinput = {};
    curAverInput = {};

    % 获取分离状态
    SPLITTED = SPLITTE_temp(i);
    xparticlei = xparticlek_1(i, :);

    % 计算裂纹尺寸参数
    a_up = xparticlei(42) - 30;
    a_down = xparticlei(22) - 30;

    % 判断裂纹阶段
    m_index = 0;
    stage = [1,2,3,4,5,6,7,0,0,0,0,8];  % 每个阶段在元胞中的位置

    while m_index == 0
        m_index = getModelIndexFunc(a_up, a_down);  % 判断裂纹阶段
        if m_index == 0
            % 随机选择其他粒子重试
            tmp_sel = randi([1, N], 1, 1);
            xparticlei = xparticlek_1(tmp_sel, :);
            a_up = xparticlei(42) - 30;
            a_down = xparticlei(22) - 30;
        end
    end

    % 根据分离状态选择代理模型
    if SPLITTED && (m_index == 3 || m_index == 5)  % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
        if m_index == 3
            curUinput = Uinput_splitted_1;
            curAverInput = averInput_splitted_1;
        elseif m_index == 5
            curUinput = Uinput_splitted_2;
            curAverInput = averInput_splitted_2;
        end
        m_name = sprintf('nn_stage%ds', stage(m_index));
    else  % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
        curUinput = Uinput_integrated{stage(m_index)};
        curAverInput = averInput_integrated{stage(m_index)};
        m_name = sprintf('nn_stage%d', stage(m_index));
    end

    % 提取粒子状态参数
    yRegSet = xparticlei(1:21);
    zRegSet = xparticlei(22:42);
    logCstar = xparticlei(43);
    gamma = xparticlei(44);

    % 调用裂纹扩展预测函数
    [yRegSet, zRegSet, SPLITTED, logCstar, gamma, deltaK_temp] = ...
        a2aNew(yRegSet, zRegSet, aver_delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step, testErrSet);

    % 存储应力强度因子增量（取各节点中的最大值）
    deltaKSet(i) = max(deltaK_temp);
end

%% ===================================================================
%% 步骤3: 计算初始时刻的失效概率
%% ===================================================================

% 使用核密度估计计算deltaK的概率密度函数
[y1, x1] = ksdensity(deltaKSet);

% 在扩展范围内重新计算密度
[y1, x1] = ksdensity(deltaKSet, linspace(min(x1), max(x1)+3, N));

% 只保留正值部分
label = find(x1 > 0);
xs = x1(label);
ys = y1(label);

% 计算断裂韧性的累积分布函数
yr = normcdf(xs, 33.4, 3.34);

% 计算积分函数
yrs = yr .* ys;

% 数值积分计算失效概率
pofCum = cumtrapz(xs, yrs);
pof = pofCum(end);

% 存储初始失效概率
PoFSet(end) = pof;

%% ===================================================================
%% 步骤4: 迭代计算直到失效概率达到临界值
%% ===================================================================

% 初始化分离标志和循环变量
SPLITTE_temp = zeros(1, N);
j = 2;
k1 = j - 1;

% 当失效概率小于临界值时继续迭代
while PoFSet(end) < 10^(-7)
    % 复制粒子状态
    xparticlek_1 = xparticleForCalII;

    % 获取当前应力范围
    aver_delta_sigma = aver_delta_sigma_set(j-1);

    % 并行计算每个粒子的裂纹扩展
    parfor i = 1:N
        % 初始化代理模型输入
        curUinput = {};
        curAverInput = {};

        % 获取分离状态
        SPLITTED = SPLITTE_temp(i);
        xparticlei = xparticlek_1(i, :);

        % 计算裂纹尺寸参数
        a_up = xparticlei(42) - 30;
        a_down = xparticlei(22) - 30;

        % 判断裂纹阶段
        m_index = 0;
        stage = [1,2,3,4,5,6,7,0,0,0,0,8];  % 每个阶段在元胞中的位置

        while m_index == 0
            m_index = getModelIndexFunc(a_up, a_down);  % 判断裂纹阶段
            if m_index == 0
                % 随机选择其他粒子重试
                tmp_sel = randi([1, N], 1, 1);
                xparticlei = xparticlek_1(tmp_sel, :);
                a_up = xparticlei(42) - 30;
                a_down = xparticlei(22) - 30;
            end
        end

        % 根据分离状态选择代理模型
        if SPLITTED && (m_index == 3 || m_index == 5)  % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
            if m_index == 3
                curUinput = Uinput_splitted_1;
                curAverInput = averInput_splitted_1;
            elseif m_index == 5
                curUinput = Uinput_splitted_2;
                curAverInput = averInput_splitted_2;
            end
            m_name = sprintf('nn_stage%ds', stage(m_index));
        else  % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
            curUinput = Uinput_integrated{stage(m_index)};
            curAverInput = averInput_integrated{stage(m_index)};
            m_name = sprintf('nn_stage%d', stage(m_index));
        end

        % 提取粒子状态参数
        yRegSet = xparticlei(1:21);
        zRegSet = xparticlei(22:42);
        logCstar = xparticlei(43);
        gamma = xparticlei(44);

        % 调用裂纹扩展预测函数
        [yRegSet, zRegSet, SPLITTED, logCstar, gamma, deltaK_temp] = ...
            a2aNew(yRegSet, zRegSet, aver_delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step, testErrSet);

        % 更新粒子状态（取实部）
        xparticleForCalII(i, :) = real([yRegSet, zRegSet, logCstar, gamma]);

        % 更新分离标志
        SPLITTE_temp(i) = SPLITTED;

        % 存储应力强度因子增量（取各节点中的最大值）
        deltaKSet(i) = max(deltaK_temp);
    end

    %% ===================================================================
    %% 计算当前时刻的失效概率
    %% ===================================================================

    % 使用核密度估计计算deltaK的概率密度函数
    [y1, x1] = ksdensity(deltaKSet);

    % 在扩展范围内重新计算密度
    [y1, x1] = ksdensity(deltaKSet, linspace(min(x1), max(x1)+3, N));

    % 只保留正值部分
    label = find(x1 > 0);
    xs = x1(label);
    ys = y1(label);

    % 计算断裂韧性的累积分布函数
    yr = normcdf(xs, 33.4, 3.34);

    % 计算积分函数
    yrs = yr .* ys;

    % 数值积分计算失效概率
    pofCum = cumtrapz(xs, yrs);
    pof = pofCum(end);

    % 添加到失效概率序列
    PoFSet = [PoFSet, pof];

    % 更新循环变量
    j = j + 1;

    % 显示计算进度
    disp(['检查间隔计算：第' num2str((j-1)*step/1950.70866) '个小时; 失效概率为' num2str(PoFSet(end))]);
end

%% ===================================================================
%% 步骤5: 计算最终检查间隔
%% ===================================================================

k2 = j - 1;
InspectionInterval = (k2 - k1) * step / 1950.70866;