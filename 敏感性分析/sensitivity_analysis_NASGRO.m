clc; clear;
close all;

%% ================== 1. 读取参数文件 ==================
filename = 'NASGRO(H-S).xlsx';   % 改成你的文件名
T = readtable(filename);
param_names = T.Properties.VariableNames;  % {'log10_D', 'p', 'K_thr', 'A'}
Xdata = table2array(T);   % N × 4

d = size(Xdata, 2);       % 参数维数 = 4

%% ================== 2. 定义参数范围 ==================
lb = min(Xdata, [], 1);
ub = max(Xdata, [], 1);

%% ================== 3. Sobol 采样（Saltelli） ==================
% rng(2026);        % 固定随机种子，保证可重复
N = 100000;         % 基样本数

A = lhsdesign(N, d);
B = lhsdesign(N, d);

A = lb + A .* (ub - lb);
B = lb + B .* (ub - lb);

% 构造 AB_i 矩阵
AB = cell(d,1);
for i = 1:d
    AB{i} = A;
    AB{i}(:,i) = B(:,i);
end

%% ================== 4. 模型输入设置 ==================
deltaK = 6;    % MPa*sqrt(m)
R = 0.5;       % 应力比
Kmax = deltaK / (1 - R);  % 最大应力强度因子

%% ================== 5. 计算模型输出 ==================
YA  = NASGRO_model(A,  deltaK, Kmax);
YB  = NASGRO_model(B,  deltaK, Kmax);

YAB = zeros(N,d);
for i = 1:d
    YAB(:,i) = NASGRO_model(AB{i}, deltaK, Kmax);
end

%% ================== 6. Sobol 指数计算 ==================
VY = var([YA; YB], 1);

S  = zeros(d,1);
ST = zeros(d,1);

for i = 1:d
    S(i)  = mean(YB .* (YAB(:,i) - YA)) / VY;
    ST(i) = 0.5 * mean((YA - YAB(:,i)).^2) / VY;
end

%% ================== 7. 输出结果 ==================
SobolTable = table(param_names', S, ST, ...
    'VariableNames', {'Parameter','FirstOrder','TotalOrder'});
disp(SobolTable)

%% ================== 8. 可视化 ==================
% Total-order
figure;
bar(ST);
set(gca,'XTick',1:d,'XTickLabel',param_names);
ylabel('Total-order Sobol index');
title('Total-order Sensitivity Index (S_T)');
grid on;

% First-order
figure;
bar(S);
set(gca,'XTick',1:d,'XTickLabel',param_names);
ylabel('First-order Sobol index');
title('First-order Sensitivity Index (S_1)');
grid on;

%% ================== 9. NASGRO 模型函数 ==================
function da_dN = NASGRO_model(X, deltaK, Kmax)
% X: N × 4 参数矩阵
% [log10(D), p, K_thr, A]

log10_D = X(:,1);
p       = X(:,2);
K_thr   = X(:,3);
A       = X(:,4);

D = 10.^log10_D;

% 避免 (deltaK - K_thr) < 0
effective_dK = max(deltaK - K_thr, 0);

% NASGRO 公式
denominator = max(1 - Kmax ./ A, 1e-6);  % 避免除零
da_dN = D .* (effective_dK).^p ./ (denominator).^(p/2);

% 防止 NaN / Inf
da_dN(~isfinite(da_dN)) = 0;

end
