%% 验证 pred_a_N 函数的测试脚本
clear; clc;

% 1. 输入从 gen_synthetic_data_spectrum.m 得到的数据
% step = 1000
t_check = [110.3872 220.7745 331.1617 441.5489 551.9362 662.3234 772.7106 883.0979 993.4851];
z = [10.9043 11.0747 12.2835 13.0287 14.2125 16.4101 17.4748 25.1182 30.2851];

% step = 10;
% t_check = [202.4079  404.8158  607.2238  809.6317 1012.0396 1214.4475 1416.8555 1619.2634 1821.6713];
% z = [10.5570 11.1956 11.9409 12.8315 13.9263 15.3394 17.3111 20.5158 30.0000];
% true_params = [-12.5479, 3.9835, 0.3326, 6.6750]; % [log_theta1, theta2, theta3, k2]

% --- 参数设置 ---
% 1. 真实参数 (用于生成参考轨迹)
theta_true = [-12.5479, 3.983501235, 0.332574645, 6.674963682];

% 2. 测试参数 (当前待验证的参数)
theta_test = [-12.3768, 3.7946, 0.4378, 6.6750];

p.theta = theta_test;
p.f = 60;                  % 应力-力转换系数 (stress_to_force)

p.k = 6.6750;      % k2
p.W = 60;                  % 试样宽度
p.B = 5;                   % 试样厚度
p.cyclesperhour = 1950.70866;

% 定义 da/dN 模型的占位符函数
% 注意：pred_a_N 传入的 delta_K 和 Kmax 已经是力对应的 K 了
% p.eq_fun = @(dk_m, theta_vec, f_factor, k_val, kmax_m) ...
%     (theta_vec(1) - theta_vec(3).*log10(theta_vec(4)) + (theta_vec(2) + theta_vec(3)).*log10(dk_m) + theta_vec(3).*(kmax_m./theta_vec(5) - 1));
p.eq_fun = @(dk_m, theta_vec, f_factor, k_val, kmax_m) ...
    (10^theta_vec(1)) * (dk_m.^theta_vec(2)) .* (max(kmax_m./theta_vec(4) - 1, 0).^theta_vec(3));
% 3. 调用 fit_a_N 进行预测
% 注意：我们需要调整 fit_a_N 内部的应力转换逻辑，或者在 fit_a_N 中加入 p.f 的使用
% 观察 fit_a_N.m 第 83-84 行：
% delta_K = (Smax - Smin) / ref_load * K_base;
% Kmax = Smax / ref_load * K_base;
% 这里 Smax 是应力，而 K_base 是基于 ref_load (力) 的。
% 所以计算出的 delta_K 物理单位不匹配。
% 我们需要修改 fit_a_N.m，使 Smax 乘以 p.f 转换为力。

fprintf('--- 开始验证 pred_a_N ---\n');

% 运行真实参数
p.theta = theta_true;
[mse_true, a_pred_true, da_hist_true] = pred_a_N(p, t_check, z);
mse_true
% 运行测试参数
p.theta = theta_test;
[mse_test, a_pred_test, da_hist_test] = pred_a_N(p, t_check, z);

% 构造增长长度矩阵: 第一行真实，第二行测试
% 注意：两者长度应一致，因为 N_check 相同
da_matrix = [da_hist_true; da_hist_test];
% 4. 可视化对比结果
figure('Color', 'w');
plot(t_check, z, 'ro', 'MarkerSize', 8, 'DisplayName', 'True Inspection (Noise Added)');
hold on;
plot(t_check, a_pred_test, 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Predicted (Test Params)');
plot(t_check, a_pred_true, 'g--', 'LineWidth', 1.2, 'DisplayName', 'Predicted (True Params)');
xlabel('Time (Flight hours)');
ylabel('Crack Length (mm)');
title(['Validation of pred\_a\_N (MSE: ', num2str(mse_test, '%.4f'), ')']);
legend('Location', 'best');
grid on;

fprintf('测试参数均方误差 MSE: %.6f\n', mse_test);

