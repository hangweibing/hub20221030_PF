function varargout = kde_sampling(mode, varargin)
%% KDE_SAMPLING 多变量核密度估计 (KDE) 参数采样主函数
%
% 参考 copula_sampling 的结构，使用 KDE 进行 4 参数分布拟合和采样。
% 特别处理：将 theta3 (列3) 和 k2 (列4) 在对数空间进行采样后再还原。
%
% 输入参数:
%   mode - 模式字符串 ('preprocess', 'sample', 'batch_sample')
%   varargin - 根据模式的可变参数
%
% 返回值根据模式而定

    persistent X_train_transformed bandwidths is_initialized
    
    % 参数索引说明:
    % 1: log_theta1_ (在原始数据中可能已经是对数形式)
    % 2: theta2
    % 3: theta3 (采用对数变换采样)
    % 4: k2     (采用对数变换采样)

    switch lower(mode)
        case 'preprocess'
            [X_train_transformed, bandwidths, is_initialized] = preprocess_kde();
            varargout = {};

        case 'sample'
            if isempty(is_initialized) || ~is_initialized
                error('KDE 未初始化，请先调用 kde_sampling(''preprocess'')');
            end
            params = sample_kde_internal(X_train_transformed, bandwidths, 1);
            varargout = {params(1), params(2), params(3), params(4)};

        case 'batch_sample'
            if isempty(is_initialized) || ~is_initialized
                error('KDE 未初始化，请先调用 kde_sampling(''preprocess'')');
            end
            n_samples = varargin{1};
            varargout{1} = sample_kde_internal(X_train_transformed, bandwidths, n_samples);

        otherwise
            error('未知的模式: %s。可用模式: preprocess, sample, batch_sample', mode);
    end
end

function [X_trans, h, is_initialized] = preprocess_kde()
%% PREPROCESS_KDE 预处理训练数据并计算 KDE 带宽

    fprintf('=== 开始 KDE 预处理 ===\n');

    % 读取训练数据
    filename = 'AM-TC4-GRO.xlsx';
    try
        T = readtable(filename);
        X_raw = table2array(T);
        [N_train, d] = size(X_raw);
        fprintf('成功读取训练数据: %d 组参数，维度 = %d\n', N_train, d);
    catch ME
        error('无法读取训练数据文件 %s: %s', filename, ME.message);
    end

    % 数据转换
    % 对 theta3 (col 3) 和 k2 (col 4) 取对数
    % 使用 max(..., eps) 防止 log(0)
    X_trans = X_raw;
    X_trans(:, 3) = log(max(X_raw(:, 3), eps));
    X_trans(:, 4) = log(max(X_raw(:, 4), eps));

    % 计算各维度的带宽 (Bandwidth)
    % 使用 MATLAB 默认的带宽选择方法
    h = zeros(1, d);
    for j = 1:d
        [~, ~, h(j)] = ksdensity(X_trans(:, j));
    end

    fprintf('KDE 带宽计算完成: [%s]\n', num2str(h, '%.4f '));

    is_initialized = true;
    fprintf('KDE 预处理完成！\n');
end

function sampled_params = sample_kde_internal(X_trans, h, n_samples)
%% SAMPLE_KDE_INTERNAL 从 KDE 分布中采样
% 实现策略：从原始转换数据中随机抽样，并叠加基于带宽的高斯扰动 (Jittering)

    [N_train, d] = size(X_trans);

    % 1. 从原始转换后的数据中随机选择索引
    idx = randi(N_train, n_samples, 1);
    base_points = X_trans(idx, :);

    % 2. 添加高斯噪声
    % 噪声标准差由带宽 h 控制
    % 简单改进：添加带宽缩减因子 (bandwidth_factor)
    % 值越小，采样点越靠近原始训练数据，分散性越低，避开不合理组合
    bandwidth_factor = 0.6; 
    noise = randn(n_samples, d) .* (h * bandwidth_factor);

    % 3. 得到变换空间的采样点
    sampled_trans = base_points + noise;

    % 4. 还原到原始物理空间
    sampled_params = sampled_trans;
    
    % 1, 2 维度保持不变 (假设 log_theta1_ 和 theta2 无需反变换)
    % 3, 4 维度通过 exp 还原 (theta3 和 k2)
    sampled_params(:, 3) = exp(sampled_trans(:, 3));
    sampled_params(:, 4) = exp(sampled_trans(:, 4));

    % 数值稳定性检查：防止极其接近 0 的负值或 NaN (虽然 exp 保证正值)
    sampled_params(:, 3) = max(sampled_params(:, 3), eps);
    sampled_params(:, 4) = max(sampled_params(:, 4), eps);
end
