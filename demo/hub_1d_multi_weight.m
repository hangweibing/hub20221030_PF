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
%% 1. 全局配置与参数设置
%% ===================================================================
% 仿真与算法参数
N = 5000;                                   % 粒子数量
step = 1000;                                % 每步循环次数
cycles_per_hour = 1950.70866;               % 每小时循环次数
R_noise = 0.5;                              % 观测噪声方差

% 参数拟合模式配置
fit_fix_mode = 'theta3';                    % 'k2': 固定 k2, 优化 theta; 'theta3': 固定 theta3, 优化 k2
k2_fixed_val = 0.394;                       % 当 fit_fix_mode='k2' 时生效
theta3_fixed_val = 0.1603;                  % 当 fit_fix_mode='theta3' 时生效

% 载荷与几何参数
spectrum_factor = 18;                       % 载荷-应力转化系数
ref_load = 100;                             % 基准载荷
W = 60;                                     % 试样宽度 (mm)
B = 5;                                      % 试样厚度 (mm)

% 真实材料参数 (Ground Truth, 用于对比)
% log_theta1=-9.7360, theta2=2.6089, theta3=0.1603, k2=0.3940
TRUE_PARAMS = [-9.7360, 2.6089, 0.1603, 0.3940]; 
true_log_theta1 = TRUE_PARAMS(1);
true_theta2 = TRUE_PARAMS(2);
true_theta3 = TRUE_PARAMS(3);
true_k2 = TRUE_PARAMS(4);

% 显示与后处理参数
display_interval = 10;                      % 信息显示间隔 (每10步显示一次)
p_up = 99.9;                                % 不确定度上边界分位数 (%)
p_low = 0.1;                                % 不确定度下边界分位数 (%)

%% ===================================================================
%% 2. 粒子滤波初始化
%% ===================================================================
n = 1;                                      % 隐性状态维度 (裂纹长度 a)

% 预采样参数组合 (log_theta1, theta2, theta3, k2)
addpath(genpath('For_PF_260131_1')) 
kde_sampling('preprocess');
fprintf('预采样 %d 个粒子的参数组合...\n', N);
particle_params = kde_sampling('batch_sample', N);
fprintf('参数预采样完成！\n');

% 可视化采样参数分布 (新增)
visualize_parameter_distribution(particle_params, 'AM-TC4-GRO.xlsx');


%% ===================================================================
%% 2. 观测数据配置 (保持原样或简化)
%% ===================================================================
% t_check = [9.1412  18.2824  27.4236  36.5648  45.7061  54.8473  63.9885  73.1297  82.2709  91.4121 100.5533 109.6945 118.8358 127.9770 137.1182 146.2594 155.4006 164.5418 173.6830];
% z = [10.7228 10.1331 10.0000 10.4430 11.6325 12.7796 13.1123 13.2576 13.8628 14.1211 15.8897 16.1198 18.7390 18.2017 19.3229 20.6861 22.0494 25.5740 30.3130];
t_check = [101.6246  203.2492  304.8738  406.4984  508.1230  609.7476  711.3722  812.9969  914.6215 1016.2461];
z = [11.0602 11.5070 12.9673 13.4202 15.3420 15.9837 18.7313 20.3204 24.5301 30.0104];

%% ===================================================================
%% 3. 数据载入与 GP 模型配置
%% ===================================================================

% 加载完整载荷谱
load('AsteixSpectraData_fake.mat', 'spectra');
spectra = repmat(reshape(spectra, 1, []), 1, 80);
spectra = spectra(2:end);

% 加载参数拟合所需的 GP 模型和人口数据
load('AM-TC4-GRO_260123_merged.mat');
load('parameter_gp_AM-TC4-GRO_combine.mat');
p_model = parameter_gp;
p_model.spectra = spectra;

% 初始化优化参数 (参考 hub20260127.m)
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

% 参数固定逻辑设置
if strcmp(fit_fix_mode, 'k2')
    p_model.fix_k2 = true; 
    p_model.k2_fixed_val = k2_fixed_val;
    p_model.fix_theta3 = false;
elseif strcmp(fit_fix_mode, 'theta3')
    p_model.fix_k2 = false; 
    p_model.fix_theta3 = true;
    p_model.theta3_fixed_val = theta3_fixed_val;
else
    p_model.fix_k2 = false;
    p_model.fix_theta3 = false;
end

% 初始化边界
theta_init = p_model.PARA0(3:end);
theta_gene = theta_init(2:end-1);
f_L_vec = zeros(p_model.numGenes, 1);
f_U_vec = zeros(p_model.numGenes, 1);
for pp = 1:p_model.numGenes
    if theta_gene(pp) > 0
        f_L_vec(pp) = 0.1; f_U_vec(pp) = 10;
    elseif theta_gene(pp) < 0
        f_L_vec(pp) = -10; f_U_vec(pp) = -0.1;
    else
        f_L_vec(pp) = 0; f_U_vec(pp) = 0;
    end
end
p_model.LB_orig = [-15; f_L_vec; 1; 0; 0];
p_model.UB_orig = [-5; f_U_vec; 10; 500; 500];

% 处理 evalstr
evalstr_tmp = p_model.evalstr1;
evalstr_tmp = regexprep(evalstr_tmp, 'c(\d+)', 'Const_pair_now($1)');
evalstr_tmp = regexprep(evalstr_tmp, 'x(\d+)', 'xtrain(:,$1)');
p_model.evalstr2 = evalstr_tmp;

% 其他必要的参数设置
p_model.parameter_K = p_model.PARA0(1:2);
p_model.theta = theta_init;
p_model.ytrain = 0.7; 
p_model.cyclesperhour = cycles_per_hour;

% 提取基础变量用于预测
f_base = p_model.theta;
k_base = p_model.parameter_K; % [delta_kth, kc]
Const_pair_now = p_model.Const_pair_now;
eq_fun = p_model.eq_fun;      % 裂纹扩展速率方程句柄

% 将关键参数封装到 p_model 以便传递给 pred_a_N
p_model.spectrum_factor = spectrum_factor;
p_model.ref_load = ref_load;
p_model.W = W;
p_model.B = B;
p_model.step_size = step;



%% ===================================================================
%% 4. 粒子群初始化 (1D 裂纹)
%% ===================================================================
% 状态向量定义: [a, log_theta1, theta2, theta3, k2] -> 暂时维持 5 维以支持 visualization, 
% 但在计算时我们将根据模型参数进行映射。
% 如果需要完全匹配 GP 模型，建议更新 kde_sampling 以采样 num_f + 2 个参数。
xparticle_curr = zeros(N, 5);

% 初始化裂纹长度
a_initial = normrnd(z(1), 0.05, N, 1);     % 以第一个观测值为基准初始化

for i = 1:N
    % 提取预采样的材料参数
    log_theta1 = particle_params(i, 1);
    theta2 = particle_params(i, 2);
    theta3 = particle_params(i, 3);
    k2 = particle_params(i, 4);
    
    % --- 物理一致性修正：强制执行参数固定逻辑 ---
    if strcmp(fit_fix_mode, 'k2')
        k2 = k2_fixed_val;       % 强制所有粒子使用固定的 k2
    elseif strcmp(fit_fix_mode, 'theta3')
        theta3 = theta3_fixed_val; % 强制所有粒子使用固定的 theta3
    end

    xparticle_curr(i, :) = [a_initial(i), log_theta1, theta2, theta3, k2];
end

% 打印状态以验证
if strcmp(fit_fix_mode, 'k2')
    fprintf('  [物理约束] 已强制固定所有粒子的 k2 = %.4f\n', k2_fixed_val);
elseif strcmp(fit_fix_mode, 'theta3')
    fprintf('  [物理约束] 已强制固定所有粒子的 theta3 = %.4f\n', theta3_fixed_val);
end

current_weight = 1/N * ones(N, 1);

%% ===================================================================
%% 5. 数据存储初始化 (预分配内存空间)
%% ===================================================================
% 计算大概需要的总步数，用于预分配空间
max_hours = t_check(end);
estimated_steps = ceil(max_hours * cycles_per_hour / step) + 100; 

Xpf = zeros(5, estimated_steps); 
a_upper = zeros(1, estimated_steps); 
a_lower = zeros(1, estimated_steps); 
POF_array = zeros(estimated_steps, 1);

mu_Kc = 33.4;
std_Kc = 3.34;

% 初始统计量 (第 1 步)
v_mask = ~isnan(xparticle_curr(:, 1));
Xpf(:, 1) = mean(xparticle_curr(v_mask, :))';
a_upper(1) = prctile(xparticle_curr(v_mask, 1), p_up);
a_lower(1) = prctile(xparticle_curr(v_mask, 1), p_low);
POF_array(1) = 0;

%% ===================================================================
%% 6. 粒子滤波主循环
%% ===================================================================
m = 2;  % 时间步索引
j = 1;  % 观测计数器
prev_PARA0_fit = []; % 用于存储拟合结果，作为下一次的起始点

total_tic = tic;
while (m-1)*step/cycles_per_hour <= t_check(end)
    iter_tic = tic;

    xparticle_prev = xparticle_curr;
    x_next_temp = zeros(N, 5);
    k_max_temp = zeros(N, 1);

    % --- 性能优化：将与粒子无关的载荷计算移出并行循环 ---
    % 1. 配置当前循环段的索引 (m=2 时起始周期为 0)
    N_now_cycle = (m-2) * step;
    idx_start = 2 * N_now_cycle + 1;
    idx_end = 2 * (N_now_cycle + step);
    
    if idx_end > length(spectra)
        error('A2A:SpectrumEnd', '载荷谱长度不足');
    end
    
    % 提取并转换载荷段 (仅执行一次，避免在大循环内重复操作)
    loads_segment = spectra(idx_start : idx_end);
    Smin_segment = loads_segment(1:2:end);
    Smax_segment = loads_segment(2:2:end);
    Pmin_segment = Smin_segment * spectrum_factor;
    Pmax_segment = Smax_segment * spectrum_factor;
    
    % 提取常量参数
    k_base_1 = k_base(1); 
    
    % 粒子预测步骤
    invalid_count = 0;
    parfor i = 1:N
        try
            % 提取当前状态
            xparticlei = xparticle_prev(i, :);

            if any(isnan(xparticlei))
                error('A2A:InvalidParticle', 'Inherited NaN');
            end

            a_curr = xparticlei(1);
            
            % 2. 计算基准应力强度因子 K_base (MPa*sqrt(mm))
            K_base = sim_K_func(a_curr, ref_load, W, B);
            
            % 3. 按比例得到 deltaK 和 Kmax (MPa*sqrt(mm))
            delta_K = (Pmax_segment - Pmin_segment) / ref_load * K_base;
            Kmax = Pmax_segment / ref_load * K_base;
            
            % 6. 转换为标准单位 (MPa*sqrt(m)) 用于公式计算
            dk_m = delta_K / sqrt(1000);
            kmax_m = Kmax / sqrt(1000);
            
            % 7. 计算 da (m/cycle) - 参考 pred_a_N.m 逻辑
            % 从粒子信息中提取当前的参数 (支持参数随观测动态更新)
            % 根据 eq_fun 公式测试：f(3) 对应 xparticlei(3), f(2) 对应 xparticlei(4)
            f_particle = [xparticlei(2); xparticlei(4); xparticlei(3)];
            k_particle = [k_base_1; xparticlei(5)];
            
            R_ratio_seg = 1 - dk_m ./ kmax_m;
            delta_K_th = k_particle(2) .* (1 - R_ratio_seg);
            idx_no_growth = (dk_m - delta_K_th) < 0;
            
            % 计算扩展速率 (log10 空间)
            da_log10 = eq_fun(dk_m, Const_pair_now, f_particle, k_particle, kmax_m);
            da_log10(idx_no_growth) = -99; % 极小值表示不扩展
            
            % 8. 累计增量并转换为 mm
            da_total_mm = sum(10.^(da_log10)) * 1e3;
            
            % 更新状态
            a_next = a_curr + da_total_mm;

            % 边界检查
            if a_next >= W || ~isreal(a_next)
                error('A2A:InvalidParticle', '无效的裂纹预测');
            end

            x_next_temp(i, :) = [a_next, xparticlei(2:end)]; % 保持参数不变
            k_max_temp(i) = kmax_m(end); % 记录最后一个周期的 Kmax 用于 POF

        catch ME
            x_next_temp(i, :) = NaN;
            k_max_temp(i) = NaN;
            invalid_count = invalid_count + 1;
        end
    end

    if invalid_count > 0 && mod(m, display_interval) == 0
        fprintf('  [预测阶段] 警告: 时间步 %d 有 %d 个粒子失效 (NaN)\n', m, invalid_count);
    end

    xparticle_curr = x_next_temp;

    % --- 记录当前步统计量 (直接写入内存数组) ---
    valid_mask = ~isnan(xparticle_curr(:, 1));
    if any(valid_mask)
        cur_Xpf = mean(xparticle_curr(valid_mask, :))';
        cur_a_up = prctile(xparticle_curr(valid_mask, 1), p_up);
        cur_a_low = prctile(xparticle_curr(valid_mask, 1), p_low);
    else
        cur_Xpf = Xpf(:, m-1);
        cur_a_up = a_upper(m-1);
        cur_a_low = a_lower(m-1);
    end
    
    % 计算失效概率 POF
    valid_K = k_max_temp(~isnan(k_max_temp));
    if ~isempty(valid_K)
        [POF_val, ~, ~, ~, stats] = calculatePOF(valid_K, mu_Kc, std_Kc);
        cur_POF = POF_val;
        
        if mod(m, display_interval) == 0
            fprintf('========== POF 计算统计信息 (步数: %d) ==========\n', m);
            fprintf('粒子数量:           %d\n', stats.num_particles);
            fprintf('K_max 均值:         %.4f MPa√m\n', stats.K_mean);
            fprintf('失效概率 POF:       %.4e\n', stats.POF);
            fprintf('======================================\n');
        end
    else
        cur_POF = POF_array(m-1);
    end
    
    Xpf(:, m) = cur_Xpf;
    a_upper(m) = cur_a_up;
    a_lower(m) = cur_a_low;
    POF_array(m) = cur_POF;

    %% 观测更新步骤
    if ((m-1)*step/cycles_per_hour < t_check(j)) && (m*step/cycles_per_hour >= t_check(j))
        fprintf('  [观测更新] 时间: %.2f h, 观测值: %.2f mm\n', t_check(j), z(j));

        % --- [Multi-Weight] 1. 参数拟合 (得到参数的"观测值") ---
        % 初始化多维度权重矩阵 [N x 5] (a, log_theta1, theta2, theta3, k2)
        multi_weights = ones(N, 5) / N; 
        
        fit_triggered = false;
        if j > 1
            fprintf('  [参数拟合] 正在进行多起始点参数拟合...\n');
            
            % 准备拟合所需的历史观测数据 (从第 1 次持续到第 j 次)
            p_model.data_a_N = [t_check(1:j)', z(1:j)'];
            p_model.num_of_data_a_N = size(p_model.data_a_N, 1);
            
            % 执行多起始点优化
            N_starts_opt = 15;
            [~, delta_kth_fit, kc_fit, theta_fit, ~, ~, pass_idx] = ...
                loss_cal_optimize_for_PF_multistart(p_model, pop_now, N_starts_opt, prev_PARA0_fit);
            
            if pass_idx == 1
                fprintf('  [参数拟合] 拟合成功！结果：f1=%.4f, f2=%.4f, f3=%.4f, kc=%.4f\n', ...
                    theta_fit(1), theta_fit(2), theta_fit(3), kc_fit);
                
                % 保存拟合结果供后续步使用
                new_fit = [delta_kth_fit, kc_fit, theta_fit'];
                prev_PARA0_fit = [prev_PARA0_fit; new_fit];
                
                % 参数观测值映射: 
                % theta_fit(1) -> log_theta1 (idx 2)
                % theta_fit(3) -> theta2     (idx 3)
                % theta_fit(2) -> theta3     (idx 4)
                % kc_fit       -> k2         (idx 5)
                param_obs = [theta_fit(1), theta_fit(3), theta_fit(2), kc_fit];
                fit_triggered = true;
                
                % 更新 p_model 初始点
                p_model.PARA0 = new_fit;
                p_model.parameter_K = p_model.PARA0(1:2);
                p_model.theta = p_model.PARA0(3:end);
            else
                fprintf('  [参数拟合] 拟合未通过约束检查，本步参数不更新权重。\n');
            end
        end

        % --- [Multi-Weight] 2. 计算各维度似然权重 ---
        % 2.1 裂纹长度权重 (a)
        a_pred = xparticle_curr(:, 1);
        for i = 1:N
            if isnan(a_pred(i))
                multi_weights(i, 1) = 1e-99;
            else
                res_a = z(j) - a_pred(i);
                multi_weights(i, 1) = (1/sqrt(2*pi*R_noise)) * exp(-0.5*(res_a^2)/R_noise) + 1e-99;
            end
        end

        % 2.2 模型参数权重 (仅在拟合成功时计算)
        if fit_triggered
            for d = 2:5
                obs_val = param_obs(d-1);
                pred_vals = xparticle_curr(:, d);
                for i = 1:N
                    if isnan(pred_vals(i))
                        multi_weights(i, d) = 1e-99;
                    else
                        res_p = obs_val - pred_vals(i);
                        multi_weights(i, d) = (1/sqrt(2*pi*R_noise)) * exp(-0.5*(res_p^2)/R_noise) + 1e-99;
                    end
                end
            end
        end

        % --- [Multi-Weight] 3. 可视化似然权重分布 (新增) ---
        if ~exist('param_obs', 'var'), param_obs = []; end
        plot_multi_weight_distribution(xparticle_curr, multi_weights, j, t_check(j), z(j), param_obs, TRUE_PARAMS);

        % --- [Multi-Weight] 4. 归一化与重采样 ---
        % 注意：重采样需在绘图之后进行，因为绘图需要展示采样前的分布状况
        xparticle_curr = resample_particles_independent(xparticle_curr, multi_weights);
        
        % 这里的 current_weight 仅作为均值计算的占位，设为 1/N
        current_weight = (1/N) * ones(N, 1);

        % 重采样后重新计算统计量并更新（确保绘图包含后验信息）
        valid_mask_res = ~isnan(xparticle_curr(:, 1));
        if any(valid_mask_res)
            Xpf(:, m) = mean(xparticle_curr(valid_mask_res, :))';
            a_upper(m) = prctile(xparticle_curr(valid_mask_res, 1), p_up);
            a_lower(m) = prctile(xparticle_curr(valid_mask_res, 1), p_low);
        end

        % --- [新增] 实时绘制裂纹长度图并保存 ---
        fig_rt = figure('Name', sprintf('Crack Length Realtime - Obs %d', j), 'Visible', 'off', 'Position', [100, 100, 1000, 600], 'Color', 'w');
        cur_time_axis = (0:m-1) * step / cycles_per_hour;
        plot(t_check(1:j), z(1:j), '^', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k'); hold on;
        plot(cur_time_axis, Xpf(1, 1:m), 'b', 'LineWidth', 3); hold on;
        plot(cur_time_axis, a_upper(1:m), 'r--', 'LineWidth', 2); hold on;
        plot(cur_time_axis, a_lower(1:m), 'r--', 'LineWidth', 2); hold off;
        xlabel('Flight hours/h', 'FontSize', 12); ylabel('Crack length/mm', 'FontSize', 12);
        legend('Experimental value', 'Prediction mean', sprintf('%.1f%% bounds', p_up - p_low), 'Location', 'best');
        title(sprintf('Crack Growth Real-time Prediction (Obs %d, Time %.2f h)', j, t_check(j)));
        grid on; set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
        
        % 确保保存文件夹存在 (使用与权重分布图一致的文件夹)
        save_path_rt = 'MultiWeight_Dist_Results';
        if ~exist(save_path_rt, 'dir'), mkdir(save_path_rt); end
        saveas(fig_rt, fullfile(save_path_rt, sprintf('CrackLength_Obs_%d.png', j)));
        close(fig_rt);

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
legend('Experimental value', 'Prediction mean', sprintf('%.1f%% bounds', p_up - p_low), 'Location', 'best');
grid on;
set(gca, 'gridlinestyle', ':', 'gridcolor', 'k');
title('1D Crack Growth Prediction with Uncertainty Bounds');

fprintf('仿真完成！总耗时: %.2f s\n', toc(total_tic));

%% ===================================================================
%% 8. 输出参数拟合历史
%% ===================================================================
if ~isempty(prev_PARA0_fit)
    fprintf('\n--- 参数拟合历史结果汇总 ---\n');
    fprintf('%-10s %-12s %-12s %-10s %-10s %-10s\n', '观测步', 'delta_Kth', 'Kc', 'f1', 'f2', 'f3');
    for row = 1:size(prev_PARA0_fit, 1)
        % 注意：j 是从 2 开始拟合的（第一次拟合是在第 2 次观测时）
        obs_idx = row + 1; 
        p_fit = prev_PARA0_fit(row, :);
        fprintf('%-10d %-12.4f %-12.4f %-10.4f %-10.4f %-10.4f\n', ...
            obs_idx, p_fit(1), p_fit(2), p_fit(3), p_fit(4), p_fit(5));
    end
    fprintf('---------------------------\n');
    
    % 输出真实参数进行对比
    % 对应关系: f1=log_theta1, f2=theta3, f3=theta2, Kc=k2
    fprintf('%-10s %-12s %-12.4f %-10.4f %-10.4f %-10.4f\n', ...
        'TRUE VAL', '---', true_k2, true_log_theta1, true_theta3, true_theta2);
    fprintf('---------------------------\n');
    
    % 可选：将历史结果导出为 Table 方便在变量浏览器查看
    para_fit_history = array2table(prev_PARA0_fit, ...
        'VariableNames', {'delta_Kth', 'Kc', 'f1', 'f2', 'f3'});
    assignin('base', 'para_fit_history', para_fit_history);
else
    fprintf('\n在此次运行中未产生成功的参数拟合结果。\n');
end
