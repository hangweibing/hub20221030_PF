%% ===================================================================
%% 1D 简化版粒子滤波裂纹扩展预测程序 (验证科研创新点)
%% 基于 hub20220916.m 简化
%% ===================================================================

% 清理工作空间
clear; close all; clc; tic; format long;

% 设置图形参数
set(groot, 'defaultFigureColor', 'white');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultAxesFontName', 'Times New Roman');

%% ===================================================================
%% 随机数种子设置
%% ===================================================================
SIM_SEED = 2023;
rng(SIM_SEED);

%% ===================================================================
%% 1. 粒子滤波算法参数配置
%% ===================================================================
n = 1;                                    % 隐性状态维度 (裂纹长度 a)
N = 5000;                                % 粒子数量

% 预采样参数组合 (log_theta1, theta2, theta3, k2)
% 注意：此处假设 kde_sampling.m 在路径中，或者已添加 main 路径
addpath('../main');
kde_sampling('preprocess');
fprintf('预采样 %d 个粒子的参数组合...\n', N);
particle_params = kde_sampling('batch_sample', N);
fprintf('参数预采样完成！\n');

%% ===================================================================
%% 采样结果可视化 (参考 hub20220916.m)
%% ===================================================================
% 加载训练数据用于对比
T = readtable('AM-TC4-GRO.xlsx');
X_train = table2array(T);
param_names = T.Properties.VariableNames;
d = size(X_train, 2);

%% 三参数联合分布可视化（三维散点图）
figure('Name', '三参数联合分布对比', 'Position', [100, 100, 1200, 800], 'Color', 'w');
param_triplets = [1,2,3; 1,2,4; 1,3,4; 2,3,4];

for k = 1:4
    idx_i = param_triplets(k,1);
    idx_j = param_triplets(k,2);
    idx_l = param_triplets(k,3);

    subplot(2,2,k);
    % 绘制训练数据的三维散点图
    scatter3(X_train(:,idx_i), X_train(:,idx_j), X_train(:,idx_l), 20, 'filled', ...
        'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
    hold on;
    % 绘制采样数据的三维散点图
    scatter3(particle_params(:,idx_i), particle_params(:,idx_j), particle_params(:,idx_l), 15, 'filled', ...
        'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');

    % 清理变量名用于显示
    clean_n_i = strrep(param_names{idx_i}, '_', '-');
    clean_n_j = strrep(param_names{idx_j}, '_', '-');
    clean_n_l = strrep(param_names{idx_l}, '_', '-');

    xlabel(clean_n_i); ylabel(clean_n_j); zlabel(clean_n_l);
    title(sprintf('%s vs %s vs %s', clean_n_i, clean_n_j, clean_n_l));
    legend('Training Data', 'Sampled', 'Location', 'best');
    grid on; view(45, 30);
end
fprintf('参数分布可视化完成！\n');

%% ===================================================================
%% 2. 观测数据配置 (保持原样或简化)
%% ===================================================================
% t_check = [9.1412  18.2824  27.4236  36.5648  45.7061  54.8473  63.9885  73.1297  82.2709  91.4121 100.5533 109.6945 118.8358 127.9770 137.1182 146.2594 155.4006 164.5418 173.6830];
% z = [10.7228 10.1331 10.0000 10.4430 11.6325 12.7796 13.1123 13.2576 13.8628 14.1211 15.8897 16.1198 18.7390 18.2017 19.3229 20.6861 22.0494 25.5740 30.3130];
t_check = [9.1412 27.4236 45.7061 63.9885 82.2709 100.5533 118.8358 137.1182 155.4006];
z = [10.7228 10.0000 11.6325 13.1123 13.8628 15.8897 18.7390 19.3229 22.0494];

R_noise = 0.5;                             % 观测噪声方差

%% ===================================================================
%% 3. 载荷谱简化 (Constant Load)
%% ===================================================================
% 直接定义载荷谱参数
Pmax = 6000;                              % 最大载荷 (N)
Pmin = 3000;                               % 最小载荷 (N)
delta_P = Pmax - Pmin;                     % 载荷范围 (N)
R_ratio = Pmin / Pmax;                     % 载荷比

% 试样几何尺寸 (用于K计算)
W = 60 ;                                    % 试样宽度 (m)
B = 5 ;                                     % 试样厚度 (m)

% 计算步长
step = 200;                               % 每步循环次数
cycles_per_hour = 1950.70866;              % 每小时循环次数

% 显示控制参数
display_interval = 100;                      % 信息显示间隔 (每100步显示一次)

%% ===================================================================
%% 4. 粒子群初始化 (1D 裂纹)
%% ===================================================================
% 状态向量定义: [a, log_theta1, theta2, theta3, k2]
xparticle_curr = zeros(N, 5);

% 初始化裂纹长度: 均值 2, 标准差 0.05
a_initial = normrnd(10, 0.05, N, 1);

for i = 1:N
    % 提取预采样的材料参数
    log_theta1 = particle_params(i, 1);
    theta2 = particle_params(i, 2);
    theta3 = particle_params(i, 3);
    k2 = particle_params(i, 4);

    xparticle_curr(i, :) = [a_initial(i), log_theta1, theta2, theta3, k2];
end

current_weight = 1/N * ones(N, 1);

%% ===================================================================
%% 5. 数据存储初始化
%% ===================================================================
% 简化版存储
Xpf = zeros(5, 150000);                  % 滤波估计值均值
a_upper = zeros(1, 150000);              % 99.9% 置信上限
a_lower = zeros(1, 150000);              % 0.1% 置信下限
POF_array = zeros(150000, 1);
mu_Kc = 33.4;
std_Kc = 3.34;

% 初始统计量
Xpf(:, 1) = mean(xparticle_curr)';
a_upper(1) = prctile(xparticle_curr(:, 1), 99);
a_lower(1) = prctile(xparticle_curr(:, 1), 1.05);

%% ===================================================================
%% 6. 粒子滤波主循环
%% ===================================================================
m = 2;  % 时间步索引
j = 1;  % 观测计数器

total_tic = tic;
while (m-1)*step/cycles_per_hour <= t_check(end)
    iter_tic = tic;

    xparticle_prev = xparticle_curr;
    x_next_temp = zeros(N, 5);
    k_max_temp = zeros(N, 1);

    % 粒子预测步骤
    invalid_count = 0;
    parfor i = 1:N
        try
            % 提取当前状态
            xparticlei = xparticle_prev(i, :);

            % 如果上一时刻已经是NaN，则直接跳过
            if any(isnan(xparticlei))
                error('A2A:InvalidParticle', 'Inherited NaN');
            end

            a_curr = xparticlei(1);
            log_theta1 = xparticlei(2);
            theta2_val = xparticlei(3);      % 避免与 loop 变量冲突
            theta3_val = xparticlei(4);
            k2_val = xparticlei(5);

            % 调用简化的 1D 扩展模型 (使用载荷而非应力)
            [da, deltaK_m] = calc_da_1d(a_curr, delta_P, R_ratio, Pmax, log_theta1, theta2_val, theta3_val, k2_val, step, W, B);

            % 更新裂纹长度
            a_next = a_curr + da;

            % 边界检查：如果裂纹长度超过宽度，标记为无效
            if a_next >= W 
                error('A2A:InvalidParticle', '裂纹长度超出试样边界(a_next=%.2f, W=%.2f)', a_next, W);
            end

            % 存储下一时刻状态
            x_next_temp(i, :) = [a_next, log_theta1, theta2_val, theta3_val, k2_val];

            % 记录 Kmax (用于 POF)
            k_max_temp(i) = deltaK_m / (1 - R_ratio);

        catch ME
            % 记录无效粒子信息
            x_next_temp(i, :) = NaN;
            k_max_temp(i) = NaN;
            invalid_count = invalid_count + 1;
        end
    end

    if invalid_count > 0 && mod(m, display_interval) == 0
        fprintf('  [预测阶段] 警告: 时间步 %d 有 %d 个粒子失效 (NaN)\n', m, invalid_count);
    end

    xparticle_curr = x_next_temp;

    % 更新均值估计
    valid_mask = ~isnan(xparticle_curr(:, 1));
    if any(valid_mask)
        Xpf(:, m) = mean(xparticle_curr(valid_mask, :))';
        a_upper(m) = prctile(xparticle_curr(valid_mask, 1), 99);
        a_lower(m) = prctile(xparticle_curr(valid_mask, 1), 1);
    else
        Xpf(:, m) = Xpf(:, m-1);
        a_upper(m) = a_upper(m-1);
        a_lower(m) = a_lower(m-1);
    end

    % 计算失效概率 POF
    valid_K = k_max_temp(~isnan(k_max_temp));
    if ~isempty(valid_K)
        [POF_val, ~, ~, ~, stats] = calculatePOF(valid_K, mu_Kc, std_Kc);
        POF_array(m) = POF_val;

        % 按指定间隔显示 POF 计算统计信息
        if mod(m, display_interval) == 0
            fprintf('========== POF 计算统计信息 (步数: %d) ==========\n', m);
            fprintf('粒子数量:           %d\n', stats.num_particles);
            fprintf('K_max 均值:         %.4f MPa√m\n', stats.K_mean);
            fprintf('K_max 标准差:       %.4f MPa√m\n', stats.K_std);
            fprintf('K_max 最大值:       %.4f MPa√m\n', stats.K_max);
            fprintf('失效概率 POF:       %.4e\n', stats.POF);
            if stats.POF > 1e-7
                fprintf('⚠️  警告: POF 超过临界值！\n');
            else
                fprintf('✓  POF 在安全范围内\n');
            end
            fprintf('======================================\n');
        end
    else
        POF_array(m) = POF_array(m-1);
    end

    %% 观测更新步骤
    if ((m-1)*step/cycles_per_hour < t_check(j)) && (m*step/cycles_per_hour >= t_check(j))
        fprintf('  [观测更新] 时间: %.2f h, 观测值: %.2f mm\n', t_check(j), z(j));

        % 计算似然权重
        a_pred = xparticle_curr(:, 1);
        weights_temp = zeros(N, 1);
        for i = 1:N
            if isnan(a_pred(i))
                weights_temp(i) = 1e-99;
            else
                res = z(j) - a_pred(i);
                weights_temp(i) = (1/sqrt(2*pi*R_noise)) * exp(-0.5*(res^2)/R_noise) + 1e-99;
            end
        end

        % 归一化
        current_weight = weights_temp ./ sum(weights_temp);

        % 更新状态估计 (仅对有效粒子进行加权平均)
        valid_idx = ~isnan(xparticle_curr(:, 1));
        if any(valid_idx)
            % 提取有效粒子的权重并重新归一化（确保权重总和为1）
            sub_weights = current_weight(valid_idx);
            normalized_sub_weights = sub_weights / sum(sub_weights);

            % 使用有效粒子更新状态估计
            Xpf(:, m) = (xparticle_curr(valid_idx, :)' * normalized_sub_weights);
        else
            % 如果极端情况下所有粒子都失效，则保持上一时刻估计
            Xpf(:, m) = Xpf(:, m-1);
        end

        % 更新观测后的置信区间 (基于加权粒子的近似处理，重采样后计算更准确)
        % 这里为了简单，先更新均值，重采样后再更新 bounds 也可以，或者在此处基于权重计算

        % 调用单独形成的重采样函数
        [xparticle_curr, current_weight] = resample_particles(xparticle_curr, current_weight);

        % 重采样后重新计算 bounds
        valid_mask_res = ~isnan(xparticle_curr(:, 1));
        if any(valid_mask_res)
            a_upper(m) = prctile(xparticle_curr(valid_mask_res, 1), 99);
            a_lower(m) = prctile(xparticle_curr(valid_mask_res, 1), 1);
        end

        j = j + 1;
    end

    m = m + 1;

    % 计算时间变量
    iter_time = toc(iter_tic);
    total_time = toc(total_tic);

    % 按指定间隔显示仿真进度信息
    if mod(m, display_interval) == 0
        disp(['已完成' num2str((m-1)*step/cycles_per_hour, '%.2f') '小时，进行了' num2str(j-1) '次观测', ...
            '，当前步耗时：' num2str(iter_time, '%.2f') 's', ...
            '，累计总耗时：' num2str(total_time, '%.2f') 's']);
    end

end

%% ===================================================================
%% 7. 结果可视化
%% ===================================================================
close all

%% 数据处理
% 提取时间轴 (m-1 个有效点)
time_axis = (0:m-2) * step / cycles_per_hour;

%% 绘制结果 - 裂纹长度图 (参考 hub20220916.m)
figure('Name', '裂纹长度预测结果', 'Position', [100, 100, 1000, 600], 'Color', 'w');
plot(t_check, z, '^', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k'); hold on;
plot(time_axis, Xpf(1, 1:m-1), 'b', 'LineWidth', 3); hold on;
plot(time_axis, a_upper(1:m-1), 'r--', 'LineWidth', 2); hold on;
plot(time_axis, a_lower(1:m-1), 'r--', 'LineWidth', 2); hold off;

% 图表格式设置
xlabel('Flight hours/h', 'FontSize', 14);
ylabel('Crack length/mm', 'FontSize', 14);
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% 图例设置
legend('Experimental value', 'Prediction mean', '99.9% bounds', 'Location', 'best');
grid on;
set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
title('1D Crack Growth Prediction with Uncertainty Bounds');

fprintf('仿真完成！总耗时: %.2f s\n', toc(total_tic));
