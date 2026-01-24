%% ===================================================================
%% 程序名称：generate_synthetic_data
%% 功能描述：生成 1D 裂纹扩展的合成观测数据（用于 PF 验证）
%% ===================================================================

% 清理工作空间
clear; close all; clc;

% 设置输入参数
a0 = 10;                           % 初始裂纹长度 (mm)
noise_std = 0.5;                   % 观测噪声标准差 (mm)

SIM_SEED = 2025;
rng(SIM_SEED);

% SIM_SEED = 2024;
% rng(SIM_SEED);

%% 1. 环境配置与参数加载
% 试样与几何参数 (同步 hub_1d_simplified.m)
W = 60;                                     % 试样宽度 (mm)
B = 5;                                      % 试样厚度 (mm)
Pmax = 6000;                                % 最大载荷 (N)
Pmin = 3000;                                 % 最小载荷 (N)
delta_P = Pmax - Pmin;                      % 载荷范围 (N)
R_ratio = Pmin / Pmax;                      % 载荷比
cycles_per_hour = 1950.70866;               % 每小时循环次数
step = 1;                                  % 仿真步长 (cycles)

% 加载真实参数库
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

fprintf('--- 合成数据生成中 ---\n');
fprintf('真实参数组合: log_theta1=%.4f, theta2=%.4f, theta3=%.4f, k2=%.4f\n', ...
    log_theta1, theta2, theta3, k2);

%% 2. 完整裂纹扩展仿真过程（用于确定总寿命）
a_curr = a0;
m = 1;
current_hour = 0;

t_full = [];
a_full_history = [];

% 设置输出间隔（每100次循环输出一次）
output_interval = 1000;

while a_curr < W / 2
    % 调用扩展速率函数
    try
        [da, ~] = calc_da_1d(a_curr, delta_P, R_ratio, Pmax, log_theta1, theta2, theta3, k2, step, W, B);
    catch
        fprintf('裂纹穿透试样或计算异常，仿真停止。\n');
        break;
    end

    % 更新状态
    a_curr = a_curr + da;
    current_hour = m * step / cycles_per_hour;

    % 记录完整轨迹
    t_full = [t_full, current_hour];
    a_full_history = [a_full_history, a_curr];

    % 每隔output_interval次循环输出当前裂纹长度
    if mod(m, output_interval) == 0
        fprintf('循环 %d: 当前裂纹长度 a = %.4f mm\n', m, a_curr);
    end

    m = m + 1;

end

% 确定总寿命
total_lifetime = t_full(end);
fprintf('总寿命: %.2f 小时\n', total_lifetime);

%% 3. 等间隔采样20个观测点
num_points = 20;
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
