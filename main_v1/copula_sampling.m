function varargout = copula_sampling(mode, varargin)
%% COPULA_SAMPLING Copula 参数采样主函数
%
% 输入参数:
%   mode - 模式字符串 ('preprocess', 'sample', 'batch_sample')
%   varargin - 根据模式的可变参数
%
% 返回值根据模式而定

    persistent X_train Rho_copula is_initialized

    switch lower(mode)
        case 'preprocess'
            [X_train, Rho_copula, is_initialized] = preprocess_copula();
            varargout = {};

        case 'sample'
            if isempty(is_initialized) || ~is_initialized
                error('Copula 未初始化，请先调用 copula_sampling(''preprocess'')');
            end
            [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = sample_single_params(X_train, Rho_copula);

        case 'batch_sample'
            if isempty(is_initialized) || ~is_initialized
                error('Copula 未初始化，请先调用 copula_sampling(''preprocess'')');
            end
            n_samples = varargin{1};
            varargout{1} = sample_batch_params(X_train, Rho_copula, n_samples);

        otherwise
            error('未知的模式: %s。可用模式: preprocess, sample, batch_sample', mode);
    end
end

function [X_train, Rho_copula, is_initialized] = preprocess_copula()
%% PREPROCESS_COPULA 预处理训练数据并拟合 Gaussian Copula

    fprintf('=== 开始 Copula 预处理 ===\n');

    % 读取训练数据
    filename = 'AM-TC4-GRO.xlsx';
    try
        T = readtable(filename);
        X_train = table2array(T);
        [N_train, d] = size(X_train);
        fprintf('成功读取训练数据: %d 组参数，维度 = %d\n', N_train, d);
    catch ME
        error('无法读取训练数据文件 %s: %s', filename, ME.message);
    end

    % Copula 预处理：边缘分布转换到 U 空间
    fprintf('进行边缘分布转换...\n');
    U_train = zeros(N_train, d);

    for j = 1:d
        try
            % 使用核密度估计得到CDF
            [f, xi] = ksdensity(X_train(:,j), 'Function','cdf');
            % 插值映射到 [0,1]
            U_train(:,j) = interp1(xi, f, X_train(:,j), 'linear', 'extrap');
        catch ME
            warning('参数 %d 的边缘分布转换失败: %s', j, ME.message);
            % 备选方案：使用经验CDF
            U_train(:,j) = tiedrank(X_train(:,j)) / (N_train + 1);
        end
    end

    % 防止数值问题
    U_train = min(max(U_train, 1e-6), 1-1e-6);

    % 拟合 Gaussian Copula
    fprintf('拟合 Gaussian Copula...\n');
    try
        Rho_copula = copulafit('Gaussian', U_train);
    catch ME
        error('Copula 拟合失败: %s', ME.message);
    end

    % 标记初始化完成
    is_initialized = true;

    fprintf('Copula 预处理完成！\n');
end

function [log_theta1_, theta2, theta3, k2] = sample_single_params(X_train, Rho_copula)
%% SAMPLE_SINGLE_PARAMS 从 Copula 空间采样单个参数组合


    % 从 Copula 空间采样
    U_sample = copularnd('Gaussian', Rho_copula, 1);

    % 验证 U 空间采样结果
    if any(isnan(U_sample)) || any(U_sample < 0) || any(U_sample > 1)
        warning('Copula 采样结果异常，使用备选方案');
        U_sample = rand(1, 4);  % 备选：均匀分布
    end

    % 映射回物理参数空间
    copula_params = zeros(1, 4);
    for j = 1:4
        try
            copula_params(j) = ksdensity(X_train(:,j), U_sample(j), 'Function','icdf');

            % 检查逆 CDF 结果是否有效
            if isnan(copula_params(j)) || isinf(copula_params(j))
                warning('参数 %d 的逆 CDF 计算失败，使用训练数据均值', j);
                copula_params(j) = mean(X_train(:,j));
            end
        catch ME
            warning('参数 %d 采样失败，使用训练数据均值: %s', j, ME.message);
            copula_params(j) = mean(X_train(:,j));
        end
    end

    % 分配到输出参数
    log_theta1_ = copula_params(1);
    theta2 = copula_params(2);
    theta3 = max(copula_params(3), min(X_train(:,3)));  % 确保 theta3 不小于训练数据最小值
    k2 = max(copula_params(4), min(X_train(:,4)));      % 确保 k2 不小于训练数据最小值


end

function params_matrix = sample_batch_params(X_train, Rho_copula, n_samples)
%% SAMPLE_BATCH_PARAMS 批量采样多个参数组合

    params_matrix = zeros(n_samples, 4);

    for i = 1:n_samples
        [params_matrix(i,1), params_matrix(i,2), params_matrix(i,3), params_matrix(i,4)] = ...
            sample_single_params(X_train, Rho_copula);
    end
end