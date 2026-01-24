function [da, deltaK_m] = calc_da_1d(a, delta_P, R, Pmax, log_theta1, theta2, theta3, k2, step, W, B)
% CALC_DA_1D Calculate 1D crack growth increment
% a:           crack length (mm)
% delta_P:     load range (N)
% R:           stress ratio
% Pmax:        maximum load (N)
% log_theta1:  model parameter (log10)
% theta2:      model parameter
% theta3:      model parameter
% k2:          model parameter
% step:        cycles per step
% W:           specimen width (mm)
% B:           specimen thickness (mm)

% 0. 边界检查
if a >= W
    error('A2A:InvalidParticle', '裂纹长度已达到或超过试样宽度 (a=%.2f, W=%.2f)', a, W);
end

% 1. 调用 SIF 计算接口
deltaK = sim_K_func(a, delta_P, W, B); % MPa*sqrt(mm)
Kmax = sim_K_func(a, Pmax, W, B);     % MPa*sqrt(mm)

% 2. 转换为标准单位 (m) 进行公式计算 (与原始 a2aFunc 逻辑保持一致)
deltaK_m = deltaK / sqrt(1000);
Kmax_m = Kmax / sqrt(1000);

% 3. 模型参数转换
theta1 = 10^(log_theta1);

% 4. 计算括号项 (Kmax/k2 - 1)^theta3
bracket_term = Kmax_m / k2 - 1;
bracket_term = max(bracket_term, 0); % 确保物理意义

% 5. 扩展速率计算 (单位: m)
da_m = step * theta1 * (deltaK_m^theta2) * (bracket_term^theta3);

% 6. 转换为 mm 输出
da = 1e3 * da_m;
end
