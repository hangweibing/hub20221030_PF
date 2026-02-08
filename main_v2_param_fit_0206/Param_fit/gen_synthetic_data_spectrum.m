%% ===================================================================
%% 程序名称：generate_synthetic_data_spectrum
%% 功能描述：基于载荷谱 AsteixSpectraData 生成 1D 裂纹扩展的合成观测数据
%% ===================================================================

% 清理工作空间
clear; close all; clc;

% 设置输入参数
a0 = 10;                           % 初始裂纹长度 (mm)
noise_std = 0.0;                   % 观测噪声标准差 (mm)
stress_to_force = 60;              % 手动设置系数：将应力转化为力 (N/MPa)
cycles_per_hour = 1950.70866;      % 每小时循环次数
step_size = 1000;                  % 计算周期 (cycles)
reference_load = 1.0;              % 基准载荷 (N)：用于计算单位载荷下的应力强度因子

SIM_SEED = 2025;
rng(SIM_SEED);

%% 1. 环境配置与参数加载
W = 60;                                     % 试样宽度 (mm)
B = 5;                                      % 试样厚度 (mm)

% 加载载荷谱 (参考 extract_sif_cycles.m#L13-14)
load('AsteixSpectraData_fake.mat', 'spectra');
spectra = repmat(reshape(spectra, 1, []), 1, 500); % 增加重复次数以防寿命较长
spectra = spectra(2:end);

% 加载真实参数库 (AM-TC4-GRO.xlsx)
T = readtable('AM-TC4-GRO.xlsx');
params_pool = table2array(T);
num_data = size(params_pool, 1);

% 随机抽取一组真实参数
idx = randi(num_data);
true_params = params_pool(idx, :);
log_theta1 = true_params(1);
theta2 = true_params(2);
theta3 = true_params(3);
k2 = true_params(4);
theta1 = 10^log_theta1;

fprintf('--- 载荷谱合成数据生成中 ---\n');
fprintf('真实参数组合: log_theta1=%.4f, theta2=%.4f, theta3=%.4f, k2=%.4f\n', ...
    log_theta1, theta2, theta3, k2);

%% 2. 完整裂纹扩展仿真过程
a_curr = a0;
N_now = 0;
current_hour = 0;

t_full = [0];
a_full_history = [a0];

fprintf('开始仿真计算...\n');

while a_curr < W / 2
    % 检查谱是否用完
    if (2 * (N_now + step_size)) > length(spectra)
        fprintf('载荷谱已用完，仿真提前结束。\n');
        break;
    end

    % 提取当前周期内的分段载荷谱
    idx_start = 2 * N_now + 1;
    idx_end = 2 * (N_now + step_size);
    loads_segment = spectra(idx_start : idx_end);
    Smin_vec = loads_segment(1:2:end);
    Smax_vec = loads_segment(2:2:end);

    % 将应力转换为力
    Pmin_vec = Smin_vec * stress_to_force;
    Pmax_vec = Smax_vec * stress_to_force;
    deltaP_vec = Pmax_vec - Pmin_vec;

    % 计算当前步的 da
    % 计算基准 K (使用参考载荷)
    K_ref_unit = sim_K_func(a_curr, reference_load, W, B);

    % 计算比例系数：K与载荷成正比
    % K_actual = K_ref_unit * (P_actual / reference_load)
    load_scaling_factor = K_ref_unit / reference_load;

    % 计算每个循环的 deltaK 和 Kmax (MPa*sqrt(mm))
    deltaK_vec = deltaP_vec * load_scaling_factor;
    Kmax_vec = Pmax_vec * load_scaling_factor;

    % 转换为标准单位 (m) 并应用模型
    dk_m = deltaK_vec / sqrt(1000);
    kmax_m = Kmax_vec / sqrt(1000);

    bracket = max(kmax_m / k2 - 1, 0);
    da_cycle_m = theta1 * (dk_m.^theta2) .* (bracket.^theta3);

    % 累加 delta_a (mm)
    delta_a_step = sum(da_cycle_m) * 1e3;

    % 更新状态
    a_curr = a_curr + delta_a_step;
    N_now = N_now + step_size;
    current_hour = N_now / cycles_per_hour;

    % 记录
    t_full = [t_full, current_hour];
    a_full_history = [a_full_history, a_curr];

    if mod(length(t_full), 50) == 0
        fprintf('当前时间: %.2f h, 裂纹长度: %.4f mm\n', current_hour, a_curr);
    end
end

total_lifetime = t_full(end);
fprintf('仿真完成。总寿命: %.2f 小时, 最终裂纹长度: %.4f mm\n', total_lifetime, a_curr);

%% 3. 等间隔采样10个观测点
num_points = 10;
t_check = linspace(0, total_lifetime, num_points);
t_check = t_check(2:end);  % 去掉起始点

% 插值得到裂纹长度
a_true_history = interp1(t_full, a_full_history, t_check, 'linear');

%% 4. 添加观测噪声
z_inspected = a_true_history + normrnd(0, noise_std, size(a_true_history));
z_inspected = max(z_inspected, a0); % 物理约束

%% 5. 可视化与输出
figure('Name', 'Synthetic Data (Spectrum)', 'Color', 'w');
plot(t_full, a_full_history, 'b-', 'LineWidth', 1.5); hold on;
plot(t_check, z_inspected, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Time (Flight hours)');
ylabel('Crack Length (mm)');
legend('True Path (Spectrum)', 'Inspections (with Noise)');
grid on;
title('Synthetic Data Generation using Load Spectrum');

fprintf('\n=== 生成的检查数据 (用于 pred_a_N 验证) ===\n');
fprintf('t_check = [%s];\n', num2str(t_check, '%.4f '));
fprintf('z = [%s];\n', num2str(z_inspected, '%.4f '));
fprintf('真实参数 (log_theta1, theta2, theta3, k2): \n[%s]\n', num2str(true_params, '%.4e '));
fprintf('应力-力转换系数: %.2f\n', stress_to_force);
fprintf('===========================================\n');

