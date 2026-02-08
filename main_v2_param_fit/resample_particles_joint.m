function xparticle_resampled = resample_particles_joint(xparticle, multi_weights, options)
% RESAMPLE_PARTICLES_JOINT 基于联合权重分布的粒子重采样
%
% 功能说明：
%   与独立重采样不同，本函数计算各维度权重的联合分布，对完整的粒子
%   （裂纹长度 + 模型参数）进行整体重采样，保持参数组合的物理一致性。
%   重采样后添加核函数扰动以维持粒子多样性。
%
% 输入参数：
%   xparticle     - [N x D] 粒子矩阵
%   multi_weights - [N x D_w] 多维度权重矩阵，每列对应一个状态维度的似然权重
%   options       - (可选) 结构体，包含以下字段：
%       .weight_mode   - 权重组合模式: 'product' (默认) 或 'weighted_sum'
%       .dim_weights   - 各维度的重要性权重 [1 x D_w]，仅在 weighted_sum 模式下使用
%       .kernel_scale  - 核函数带宽缩放因子 (默认 0.1)
%
% 输出参数：
%   xparticle_resampled - [N x D] 重采样后的粒子矩阵

%% 参数解析
[N, D] = size(xparticle);
[~, Dw] = size(multi_weights);

% 默认选项
if nargin < 3 || isempty(options)
    options = struct();
end

weight_mode = get_option(options, 'weight_mode', 'product');
dim_weights = get_option(options, 'dim_weights', ones(1, Dw));
kernel_scale = get_option(options, 'kernel_scale', 0.1);

%% 1. 计算联合权重
switch weight_mode
    case 'product'
        % 乘积模式：各维度权重相乘（对数空间求和以提高数值稳定性）
        log_weights = sum(log(multi_weights + 1e-300), 2); % 避免 log(0)
        joint_weights = exp(log_weights - max(log_weights)); % 数值稳定化
        
    case 'weighted_sum'
        % 加权求和模式：各维度权重加权平均
        joint_weights = multi_weights * dim_weights(:);
        
    otherwise
        error('不支持的权重模式: %s', weight_mode);
end

% 归一化联合权重
joint_weights = joint_weights / sum(joint_weights);

%% 2. 计算有效粒子数 (Effective Sample Size)
Neff = 1 / sum(joint_weights.^2);
fprintf('  [联合重采样] 有效粒子数: %.1f / %d (%.1f%%)\n', Neff, N, 100*Neff/N);

%% 3. 系统重采样 (Systematic Resampling)
% 使用系统重采样算法，相比多项式重采样方差更小
cumsum_weights = cumsum(joint_weights);
u0 = rand() / N; % 初始随机偏移
u = u0 + (0:N-1)' / N; % 均匀间隔的采样点

% 重采样索引
resample_idx = zeros(N, 1);
j = 1;
for i = 1:N
    while u(i) > cumsum_weights(j)
        j = j + 1;
    end
    resample_idx(i) = j;
end

xparticle_resampled = xparticle(resample_idx, :);

%% 5. 核函数扰动 (Kernel Density Perturbation)
% 由于程序处理的是二维裂纹，对裂纹节点坐标随意添加扰动容易导致裂纹面超出物理边界
% 因此仅对模型参数（维度 43-46）施加高斯扰动，不对裂纹几何坐标（维度 1-42）施加扰动

valid_mask = ~isnan(xparticle_resampled);
bandwidth = zeros(1, D);

% 确定需要施加扰动的维度范围
if D == 46
    % 46维情况：仅对参数维度 (43:46) 计算带宽
    param_dims = 43:46;
    fprintf('  [核扰动] 检测到46维粒子，仅对参数维度 (43-46) 施加扰动\n');
elseif D == 5
    % 5维情况（1D裂纹）：对除裂纹长度外的参数维度 (2:5) 施加扰动
    param_dims = 2:5;
    fprintf('  [核扰动] 检测到5维粒子，对参数维度 (2-5) 施加扰动\n');
else
    % 其他情况：对所有维度施加扰动（保守策略）
    param_dims = 1:D;
    fprintf('  [核扰动] 未知维度配置，对所有维度施加扰动\n');
end

% 使用 Silverman's rule of thumb 估计带宽（仅针对参数维度）
for d = param_dims
    valid_data = xparticle_resampled(valid_mask(:, d), d);
    if ~isempty(valid_data)
        sigma_d = std(valid_data);
        n_valid = length(valid_data);
        % Silverman's rule: h = sigma * (4/(3*n))^(1/5) ≈ 1.06 * sigma * n^(-1/5)
        bandwidth(d) = 1.06 * sigma_d * n_valid^(-0.2) * kernel_scale;
    else
        bandwidth(d) = 0;
    end
end

% 添加高斯核扰动（仅对参数维度）
for d = param_dims
    if bandwidth(d) > 0
        perturbation = randn(N, 1) * bandwidth(d);
        xparticle_resampled(:, d) = xparticle_resampled(:, d) + perturbation;
    end
end

% 输出带宽信息
if D == 46
    fprintf('  [联合重采样] 完成！参数维度带宽: [%.4e, %.4e, %.4e, %.4e]\n', ...
        bandwidth(43), bandwidth(44), bandwidth(45), bandwidth(46));
elseif D == 5
    fprintf('  [联合重采样] 完成！参数维度带宽: [%.4e, %.4e, %.4e, %.4e]\n', ...
        bandwidth(2), bandwidth(3), bandwidth(4), bandwidth(5));
else
    fprintf('  [联合重采样] 完成！\n');
end

end

%% 辅助函数：获取选项值
function val = get_option(options, field, default)
    if isfield(options, field)
        val = options.(field);
    else
        val = default;
    end
end
