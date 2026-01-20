%% 测试脚本：比较多元高斯分布和Gaussian Copula两种采样方法
clc; clear;

% 设置警告和图形参数
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
set(groot, 'defaultAxesTickLabelInterpreter', 'none');
set(groot, 'defaultLegendInterpreter', 'none');

%% 1. 读取数据
filename = 'AM-TC4-GRO.xlsx';
T = readtable(filename);

% 转为数值矩阵 (N x 4)
X = table2array(T);
param_names = T.Properties.VariableNames;

[N, d] = size(X);
fprintf('读取到 %d 组参数，维度 = %d\n', N, d);

%% 2. 估计均值和协方差（用于多元高斯）
mu = mean(X, 1);          % 1 x d
Sigma = cov(X);           % d x d

% 检查协方差矩阵是否正定
[eigvec, eigval] = eig(Sigma);
if any(diag(eigval) <= 0)
    warning('协方差矩阵非正定，可能需要正则化');
end

%% =============================================================================
%% 方法1: 多元高斯联合分布采样
%% =============================================================================
fprintf('\n=== 方法1: 多元高斯联合分布采样 ===\n');

% 抽样
Ns = 10000;   % 抽样数量
X_sample_gauss = mvnrnd(mu, Sigma, Ns);

% 可视化检查
figure('Name', '多元高斯采样结果');
for i = 1:d
    subplot(2,2,i)
    histogram(X_sample_gauss(:,i), 50, 'Normalization','pdf');
    hold on
    histogram(X(:,i), 15, 'Normalization','pdf');
    clean_name = strrep(param_names{i}, '_', '-');
    title(sprintf('%s - 高斯采样', clean_name), 'Interpreter', 'none')
    legend('Sampled','Original')
end

fprintf('多元高斯采样完成，生成 %d 个样本\n', Ns);

%% =============================================================================
%% 方法2: Gaussian Copula 联合分布采样
%% =============================================================================
fprintf('\n=== 方法2: Gaussian Copula 联合分布采样 ===\n');

%% 2.1 边缘分布：经验CDF (KDE)
U = zeros(N, d);

for j = 1:d
    % 使用核密度估计得到CDF
    [f, xi] = ksdensity(X(:,j), 'Function','cdf');
    % 插值映射到 [0,1]
    U(:,j) = interp1(xi, f, X(:,j), 'linear', 'extrap');
end

% 防止数值问题
U = min(max(U, 1e-6), 1-1e-6);

%% 2.2 拟合 Gaussian Copula 的相关矩阵
Rho = copulafit('Gaussian', U);
fprintf('Copula 相关矩阵：\n');
disp(Rho);

%% 2.3 Copula 空间采样
U_sample = copularnd('Gaussian', Rho, Ns);

%% 2.4 从 U 空间映射回物理参数空间
X_sample_copula = zeros(Ns, d);

for j = 1:d
    % 逆CDF（仍然用 KDE）
    X_sample_copula(:,j) = ksdensity(X(:,j), U_sample(:,j), ...
                              'Function','icdf');
end

%% 2.5 可视化对比
figure('Name', 'Gaussian Copula采样结果');
for i = 1:d
    subplot(2,2,i)
    histogram(X_sample_copula(:,i), 50, 'Normalization','pdf');
    hold on
    histogram(X(:,i), 15, 'Normalization','pdf');
    clean_name = strrep(param_names{i}, '_', '-');
    title(sprintf('%s - Copula采样', clean_name), 'Interpreter', 'none')
    legend('Copula Sampled','Original')
end

fprintf('Gaussian Copula 采样完成，生成 %d 个样本\n', Ns);

%% =============================================================================
%% 方法对比分析
%% =============================================================================
fprintf('\n=== 方法对比分析 ===\n');

% 计算样本统计量
fprintf('原始数据统计量:\n');
for i = 1:d
    fprintf('%s: 均值=%.4f, 方差=%.4f\n', param_names{i}, mean(X(:,i)), var(X(:,i)));
end

fprintf('\n多元高斯采样统计量:\n');
for i = 1:d
    fprintf('%s: 均值=%.4f, 方差=%.4f\n', param_names{i}, mean(X_sample_gauss(:,i)), var(X_sample_gauss(:,i)));
end

fprintf('\nGaussian Copula采样统计量:\n');
for i = 1:d
    fprintf('%s: 均值=%.4f, 方差=%.4f\n', param_names{i}, mean(X_sample_copula(:,i)), var(X_sample_copula(:,i)));
end

% 计算KS检验（分布相似性检验）
fprintf('\nKS检验结果 (p值越接近1，分布越相似):\n');
for i = 1:d
    [~, p_gauss] = kstest2(X(:,i), X_sample_gauss(:,i));
    [~, p_copula] = kstest2(X(:,i), X_sample_copula(:,i));
    fprintf('%s: 高斯采样 p=%.4f, Copula采样 p=%.4f\n', param_names{i}, p_gauss, p_copula);
end

%% =============================================================================
%% 联合分布可视化（二维投影）
%% =============================================================================
figure('Name', '联合分布对比');

% 计算需要显示的参数对数量 (4个参数的两两组合 = 6对)
param_pairs = [
    1,2; 1,3; 1,4;  % 第一行：参数1与其他参数的组合
    2,3; 2,4; 3,4   % 第二行：其余参数的组合
];

for k = 1:6
    i = param_pairs(k,1);
    j = param_pairs(k,2);

    subplot(2,3,k);
    scatter(X(:,i), X(:,j), 10, 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    scatter(X_sample_gauss(:,i), X_sample_gauss(:,j), 5, 'r', 'filled', 'MarkerFaceAlpha', 0.2);
    scatter(X_sample_copula(:,i), X_sample_copula(:,j), 5, 'g', 'filled', 'MarkerFaceAlpha', 0.2);

    % 清理变量名用于显示（移除下划线等特殊字符）
    clean_name_i = strrep(param_names{i}, '_', '-');
    clean_name_j = strrep(param_names{j}, '_', '-');

    xlabel(clean_name_i, 'Interpreter', 'none');
    ylabel(clean_name_j, 'Interpreter', 'none');
    legend('原始', '高斯', 'Copula');
    title(sprintf('%s vs %s', clean_name_i, clean_name_j), 'Interpreter', 'none');
end

fprintf('\n测试脚本执行完成！\n');
fprintf('生成的三张图：\n');
fprintf('1. 多元高斯采样结果\n');
fprintf('2. Gaussian Copula采样结果\n');
fprintf('3. 联合分布对比\n');