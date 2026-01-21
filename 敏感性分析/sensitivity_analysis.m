clc; clear;
close all;

%% ================== 1. 读取参数文件 ==================
filename = 'AM-TC4-GRO.xlsx';
T = readtable(filename);
param_names = T.Properties.VariableNames;
Xdata = table2array(T);   % N × 4

d = size(Xdata, 2);       % 参数维数 = 4

%% ================== 2. 定义参数范围 ==================
lb = min(Xdata, [], 1);
ub = max(Xdata, [], 1);

%% ================== 3. 基于 Sobol 序列的高级采样 (优化点) ==================
N = 10000; % 基样本数 (Sobol序列对于此类模型通常 10000 即可获得高精度)

% 创建 2*d 维度的 Sobol 序列对象
% 使用 'Skip' 和 'Leap' 避开序列初期的相关性
p = sobolset(2*d, 'Skip', 1000, 'Leap', 100);
p = scramble(p, 'MatousekAffineOwen'); % 引入打乱机制增强随机分布性
R = net(p, N); % 提取 N 个样本点 [N x 2d]

% 将样本点缩放到实际参数范围 [lb, ub]
VarMin = [lb, lb]; 
VarMax = [ub, ub];
R_scaled = VarMin + R .* (VarMax - VarMin);

% 拆分为 A 矩阵和 B 矩阵
A = R_scaled(:, 1:d);
B = R_scaled(:, d+1:end);

%% ================== 4. 模型输入设置 ==================
deltaK = 6;    % MPa*sqrt(m)
R_ratio = 0.8; % 应力比

%% ================== 5. 计算模型输出 ==================
YA = crack_model_PSR(A, deltaK, R_ratio);
YB = crack_model_PSR(B, deltaK, R_ratio);

% 预分配敏感性指数
S  = zeros(d, 1);
ST = zeros(d, 1);
VY = var([YA; YB], 1); % 方差分母

for i = 1:d
    % 构造 AB_i 矩阵：取 A，但将第 i 列替换为 B 的第 i 列
    AB_i = A;
    AB_i(:, i) = B(:, i);
    
    % 执行模型
    YAB_i = crack_model_PSR(AB_i, deltaK, R_ratio);
    
    % Sobol 指数估算 (使用 Saltelli 和 Jansen 建议的高稳定性算子)
    S(i)  = mean(YB .* (YAB_i - YA)) / VY;            % 一阶敏感性 (参数独立贡献)
    ST(i) = 0.5 * mean((YA - YAB_i).^2) / VY;        % 全阶敏感性 (含耦合贡献)
end

%% ================== 6. 输出结果 ==================
SobolTable = table(param_names', S, ST, ...
    'VariableNames', {'Parameter','FirstOrder','TotalOrder'});
disp('--- Sobol 敏感性分析结果 ---')
disp(SobolTable)

%% ================== 7. 可视化优化 ==================
figure('Color', 'w', 'Name', 'Crack Model Sensitivity Analysis');

% 子图 1：一阶敏感性
subplot(1, 2, 1);
bar(S, 'FaceColor', [0.3 0.6 0.9]);
set(gca, 'XTick', 1:d, 'XTickLabel', param_names);
xlabel('Parameters');
ylabel('First-order Index (S_1)');
title('First-order Sensitivity');
grid on;

% 子图 2：全阶敏感性
subplot(1, 2, 2);
bar(ST, 'FaceColor', [0.9 0.4 0.4]);
set(gca, 'XTick', 1:d, 'XTickLabel', param_names);
xlabel('Parameters');
ylabel('Total-order Index (S_T)');
title('Total-order Sensitivity');
grid on;

% 整体美化
sgtitle(['Sobol Sensitivity Analysis (N = ', num2str(N), ')']);
