%% ===================================================================
%% 测试脚本：calculatePOF 函数测试
%% 功能描述：测试基于粒子滤波的失效概率计算功能
%% ===================================================================

clear all;
clc;

fprintf('========== 开始测试 calculatePOF 函数 ==========\n');

%% ===================================================================
%% 步骤 1: 准备测试数据
%% ===================================================================

% 模拟粒子滤波结果：各粒子的最大应力强度因子 (MPa√m)
% 这里假设有100个粒子，K值服从正态分布
num_particles = 100;
mean_K = 15.0;      % 应力强度因子的均值
std_K = 3.0;        % 应力强度因子的标准差

% 生成随机应力强度因子数据
particles_K_max = mean_K + std_K * randn(num_particles, 1);

% 断裂韧性参数 (典型钢材参数)
mu_Kc = 33.4;       % 断裂韧性均值 (MPa√m)
std_Kc = 3.34;      % 断裂韧性标准差 (MPa√m)

fprintf('测试数据准备完成：\n');
fprintf('  粒子数量: %d\n', num_particles);
fprintf('  K_max 均值: %.2f MPa√m\n', mean(particles_K_max));
fprintf('  K_max 标准差: %.2f MPa√m\n', std(particles_K_max));
fprintf('  K_c 均值: %.2f MPa√m\n', mu_Kc);
fprintf('  K_c 标准差: %.2f MPa√m\n', std_Kc);
fprintf('\n');

%% ===================================================================
%% 步骤 2: 调用 calculatePOF 函数
%% ===================================================================

fprintf('正在计算失效概率...\n');

try
    % 基本用法：只传入必需参数
    [POF, f_K_vals, K_grid, F_Kc_vals] = calculatePOF(particles_K_max, mu_Kc, std_Kc);

    fprintf('✓ 计算成功！\n\n');

    %% ===================================================================
    %% 步骤 3: 显示结果
    %% ===================================================================

    fprintf('========== 计算结果 ==========\n');
    fprintf('失效概率 POF: %.6f (%.2e)\n', POF, POF);
    fprintf('临界 POF (10^-7): %.2e\n', 1e-7);

    % 判断是否超过临界值
    if POF > 1e-7
        fprintf('⚠️  警告: POF 超过临界值！结构可能存在失效风险。\n');
    else
        fprintf('✓  POF 在安全范围内，结构相对安全。\n');
    end

    fprintf('\n');

    % 显示统计信息
    fprintf('========== 统计信息 ==========\n');
    fprintf('网格点数量: %d\n', length(K_grid));
    fprintf('K 范围: [%.2f, %.2f] MPa√m\n', min(K_grid), max(K_grid));
    fprintf('f_K 峰值: %.4f\n', max(f_K_vals));
    fprintf('F_Kc 范围: [%.4f, %.4f]\n', min(F_Kc_vals), max(F_Kc_vals));

catch ME
    fprintf('✗ 计算失败！\n');
    fprintf('错误信息: %s\n', ME.message);
    fprintf('错误位置: %s (第%d行)\n', ME.stack(1).file, ME.stack(1).line);
end

fprintf('\n========== 测试完成 ==========\n');

%% ===================================================================
%% 可选：绘制结果图（取消注释以启用）
%% ===================================================================

% figure;
% subplot(2,2,1);
% histogram(particles_K_max, 20, 'Normalization', 'pdf');
% hold on;
% plot(K_grid, f_K_vals, 'r-', 'LineWidth', 2);
% title('应力强度因子分布');
% xlabel('K (MPa√m)');
% ylabel('概率密度');
% legend('直方图', 'KDE拟合');
% grid on;
%
% subplot(2,2,2);
% plot(K_grid, F_Kc_vals, 'b-', 'LineWidth', 2);
% title('断裂韧性累积分布');
% xlabel('K (MPa√m)');
% ylabel('累积概率');
% grid on;
%
% subplot(2,2,3);
% integrand = f_K_vals .* F_Kc_vals;
% plot(K_grid, integrand, 'g-', 'LineWidth', 2);
% title('积分函数 f_K(K) × F_{Kc}(K)');
% xlabel('K (MPa√m)');
% ylabel('积分值');
% grid on;
%
% subplot(2,2,4);
% text(0.1, 0.8, sprintf('失效概率 POF = %.2e', POF), 'FontSize', 12);
% text(0.1, 0.6, sprintf('粒子数 = %d', num_particles), 'FontSize', 10);
% text(0.1, 0.4, sprintf('K_{max} 均值 = %.2f', mean(particles_K_max)), 'FontSize', 10);
% text(0.1, 0.2, sprintf('K_c 均值 = %.2f', mu_Kc), 'FontSize', 10);
% axis off;
% title('计算结果汇总');