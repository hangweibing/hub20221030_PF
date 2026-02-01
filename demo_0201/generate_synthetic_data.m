%% ===================================================================
%% 程序名称：generate_synthetic_data
%% 功能描述：生成 1D 裂纹扩展的合成观测数据（用于 PF 验证）
%% ===================================================================

% 清理工作空间
clear; close all; clc;

% 设置输入参数
a0 = 10;                           % 初始裂纹长度 (mm)
noise_std = 0.0;                   % 观测噪声标准差 (mm)

SIM_SEED = 2025;
rng(SIM_SEED);

% SIM_SEED = 2024;
% rng(SIM_SEED);

%% 1. 环境配置与参数加载
% 路径配置
addpath(genpath('For_PF_260127_1'));
% addpath('../main');

% 试样与几何参数 (同步 hub_1d_simplified.m)
W = 60;                                     % 试样宽度 (mm)
B = 5;                                      % 试样厚度 (mm)
spectrum_factor = 18;                       % 载荷-应力转化系数 (Stress to Force)
ref_load = 100;                             % 基准载荷 (N)
cycles_per_hour = 1950.70866;               % 每小时循环次数
step = 100;                                 % 仿真步长 (cycles)

% 加载完整载荷谱
load('AsteixSpectraData_fake.mat', 'spectra');
spectra = repmat(reshape(spectra, 1, []), 1, 80);
spectra = spectra(2:end);

% 加载参数拟合所需的 GP 模型和人口数据
load('AM-TC4-GRO_260123_merged.mat');
load('parameter_gp_AM-TC4-GRO_combine.mat');
p_model = parameter_gp;
PARA0 = pop_now{1,6};
k_base = PARA0(1:2); % [delta_kth, kc]
Const_pair_now = p_model.Const_pair_now;
eq_fun = p_model.eq_fun;      % 裂纹扩展速率方程句柄

% 加载真实参数库
T = readtable('AM-TC4-GRO.xlsx');
params_pool = table2array(T);
num_data = size(params_pool, 1);

% 手动选取真实参数行 (注意：1 是标题，2 是第一组数据)
target_excel_row = 27; 

true_params = params_pool(target_excel_row, :);
log_theta1 = true_params(1);
theta2 = true_params(2);
theta3 = true_params(3);
k2_val = true_params(4);

fprintf('--- 合成数据生成中 ---\n');
fprintf('真实参数组合: log_theta1=%.4f, theta2=%.4f, theta3=%.4f, k2=%.4f\n', ...
    log_theta1, theta2, theta3, k2_val);

%% 2. 完整裂纹扩展仿真过程（用于确定总寿命）
a_curr = a0;
m = 1; % 仿真步数
N_total_cycles = 0; % 累计循环数

t_full = [0];
a_full_history = [a0];

% 设置输出间隔
output_interval = 10;

fprintf('开始仿真扩展...\n');
while a_curr < W / 2
    % 1. 配置当前循环段的索引
    N_start_cycle = N_total_cycles;
    
    % 2. 计算基准应力强度因子 K_base (MPa*sqrt(mm))
    K_base = sim_K_func(a_curr, ref_load, W, B);
    
    % 3. 提取载荷谱循环
    idx_start = 2 * N_start_cycle + 1;
    idx_end = 2 * (N_start_cycle + step);
    
    if idx_end > length(spectra)
        fprintf('载荷谱长度不足，仿真停止。\n');
        break;
    end
    
    loads_segment = spectra(idx_start : idx_end);
    Smin_segment = loads_segment(1:2:end);
    Smax_segment = loads_segment(2:2:end);
    
    % 4. 转化为外载荷 (Force)
    Pmin_segment = Smin_segment * spectrum_factor;
    Pmax_segment = Smax_segment * spectrum_factor;
    
    % 5. 按比例得到 deltaK 和 Kmax (MPa*sqrt(mm))
    delta_K = (Pmax_segment - Pmin_segment) / ref_load * K_base;
    Kmax = Pmax_segment / ref_load * K_base;
    
    % 6. 转换为标准单位 (MPa*sqrt(m)) 用于公式计算
    dk_m = delta_K / sqrt(1000);
    kmax_m = Kmax / sqrt(1000);
    
    % 7. 计算 da (m/cycle)
    % 映射关系: f = [log_theta1; theta3; theta2], k = [delta_kth; k2]
    f_particle = [log_theta1; theta3; theta2];
    k_particle = [k_base(1); k2_val];
    
    R_ratio_seg = 1 - dk_m ./ kmax_m;
    delta_K_th_seg = k_particle(2) .* (1 - R_ratio_seg);
    idx_no_growth = (dk_m - delta_K_th_seg) < 0;
    
    % 计算扩展速率 (log10 空间)
    da_log10 = eq_fun(dk_m, Const_pair_now, f_particle, k_particle, kmax_m);
    da_log10(idx_no_growth) = -99; % 极小值表示不扩展
    
    % 8. 累计增量并转换为 mm
    da_total_mm = sum(10.^(da_log10)) * 1e3;
    
    % 更新状态
    a_curr = a_curr + da_total_mm;
    N_total_cycles = N_total_cycles + step;
    current_hour = N_total_cycles / cycles_per_hour;

    % 记录完整轨迹
    t_full = [t_full, current_hour];
    a_full_history = [a_full_history, a_curr];

    % 每隔 output_interval 步输出当前裂纹长度
    if mod(m, output_interval) == 0
        fprintf('步数 %d (循环 %d): 当前裂纹长度 a = %.4f mm\n', m, N_total_cycles, a_curr);
    end

    m = m + 1;
    
    % 物理边界检查
    if a_curr >= W
        break;
    end
end

% 确定总寿命
total_lifetime = t_full(end);
fprintf('总寿命: %.2f 小时\n', total_lifetime);

%% 3. 等间隔采样观测点
num_points = 11;
t_check = linspace(0, total_lifetime, num_points);
t_check = t_check(2:end);  % 去掉t=0的点

% 对每个观测时间点进行插值得到裂纹长度
a_true_history = interp1(t_full, a_full_history, t_check, 'linear');

%% 4. 添加观测噪声
% z = a_true + N(0, noise_std^2)
z_inspected = a_true_history + normrnd(0, noise_std, size(a_true_history));

% 确保 z 不小于 a0 (物理常识)
z_inspected = max(z_inspected, a0);

%% 5. 可视化生成的数据
figure('Name', 'Synthetic Data Generation', 'Color', 'w');
plot(t_full, a_full_history, 'b-', 'LineWidth', 1); hold on;
plot(t_check, z_inspected, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('Time (Flight hours)');
ylabel('Crack Length (mm)');
legend('True Growth Path', 'Synthetic Inspections (with Noise)');
grid on;
title('Synthetic Data for 1D Crack Growth');

fprintf('生成的观测点数量: %d\n', length(t_check));
fprintf('最后观测时间: %.2f 小时\n', t_check(end));
fprintf('----------------------\n');

%% 6. 输出结果（hub_1d_simplified.m 格式）
fprintf('=== 合成数据生成完成 ===\n');
% 观测时间点 (小时)
fprintf('t_check = [%s];\n', num2str(t_check, '%.4f '));
% 观测值 (裂纹长度)
fprintf('z = [%s];\n', num2str(z_inspected, '%.4f '));
fprintf('真实参数: [%s]\n', num2str(true_params, '%.4f '));
fprintf('========================\n');
