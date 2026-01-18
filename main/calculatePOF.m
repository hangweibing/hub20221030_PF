%% ===================================================================
%% 函数名称：calculatePOF
%% 功能描述：基于粒子滤波的失效概率计算
%% ===================================================================
function [POF, f_K_vals, K_grid, F_Kc_vals] = calculatePOF(particles_K_max, mu_Kc, std_Kc, varargin)
% CALCULATEPOF 根据应力-强度干涉理论计算失效概率
%
% 理论依据：
%   当应力强度因子 K 超过断裂韧性 K_c 时，结构发生失效
%   POF = ∫_{-∞}^{∞} f_K(K) [∫_{-∞}^K f_{K_c}(K_c) dK_c] dK
%      = ∫_{-∞}^{∞} f_K(K) · F_{K_c}(K) dK
%
% 输入参数：
%   particles_K_max: 各粒子的最大应力强度因子向量 [N×1]
%   mu_Kc:          断裂韧性的均值 (MPa√m)
%   std_Kc:         断裂韧性的标准差 (MPa√m)
%   varargin:       可选参数
%                   - 'NumGridPoints': KDE网格点数量，默认100
%                   - 'Kernel':        核函数类型，默认'normal'
%                   - 'KdensitySupport': 支撑域类型，默认'unbounded'
%
% 输出参数：
%   POF:            失效概率 (0~1之间)
%   f_K_vals:       拟合的 f_K(K) 在网格点上的值
%   K_grid:         KDE的网格点
%   F_Kc_vals:      F_{K_c}(K) 在网格点上的值
%
%% ===================================================================

%% 输入验证
if nargin < 3
    error('calculatePOF: 至少需要3个输入参数 (particles_K_max, mu_Kc, std_Kc)');
end

% 确保 particles_K_max 是列向量
if isrow(particles_K_max)
    particles_K_max = particles_K_max';
end

% 移除 NaN、Inf 值，并只保留大于0的K值
valid_idx = isfinite(particles_K_max) & (particles_K_max > 0);
particles_K_max = particles_K_max(valid_idx);


%% 可选参数解析
p = inputParser;
addParameter(p, 'NumGridPoints', 100, @isnumeric);      % KDE网格点数
addParameter(p, 'Kernel', 'normal', @ischar);            % 核函数类型
addParameter(p, 'KdensitySupport', 'unbounded', @ischar); % 支撑域类型
parse(p, varargin{:});

num_grid_points = p.Results.NumGridPoints;
kernel_type = p.Results.Kernel;
kde_support = p.Results.KdensitySupport;

%% ===================================================================
%% 步骤 1: 核密度估计 (KDE) 拟合 f_K(K)
%% ===================================================================

% 计算 KDE（自动选择最优带宽）
% 'ksdensity' 返回概率密度函数在指定网格点上的值
try
    [f_K_vals, K_grid] = ksdensity(particles_K_max, ...
                                    'NumPoints', num_grid_points, ...
                                    'Kernel', kernel_type, ...
                                    'Support', kde_support);
catch ME
    warning(ME.identifier, 'calculatePOF: KDE计算失败 (%s)，使用默认参数重试', ME.message);
    [f_K_vals, K_grid] = ksdensity(particles_K_max, 'NumPoints', num_grid_points);
end

% 确保 f_K 为非负（理论上应该满足）
f_K_vals(f_K_vals < 0) = 0;

%% ===================================================================
%% 步骤 2: 计算断裂韧性 K_c 的累积分布函数 F_{K_c}(K)
%% ===================================================================

% 对于正态分布 K_c ~ N(mu_Kc, std_Kc^2)，使用 normcdf
F_Kc_vals = normcdf(K_grid, mu_Kc, std_Kc);

%% ===================================================================
%% 步骤 3: 计算失效概率 POF (应力-强度干涉理论)
%% ===================================================================

% 使用梯形积分法计算：
% POF = ∫ f_K(K) · F_{K_c}(K) dK
integrand = f_K_vals .* F_Kc_vals;

% 数值积分
POF = trapz(K_grid, integrand);

% 确保 POF 在 [0, 1] 范围内（理论上应该满足）
POF = max(0, min(1, POF));

%% ===================================================================
%% 诊断信息（调试用）
%% ===================================================================

% 如果 POF 异常，输出警告
if POF > 0.5
    warning('calculatePOF: POF = %.4f 异常偏高，请检查输入数据', POF);
end

% 可选：输出统计信息（取消注释以启用）
fprintf('========== POF 计算统计信息 ==========\n');
% fprintf('粒子数量:           %d\n', length(particles_K_max));
fprintf('K_max 均值:         %.4f MPa√m\n', mean(particles_K_max) / sqrt(1000));
fprintf('K_max 标准差:       %.4f MPa√m\n', std(particles_K_max) / sqrt(1000));
fprintf('K_max 最大值:       %.4f MPa√m\n', max(particles_K_max) / sqrt(1000));
fprintf('失效概率 POF:       %.4e\n', POF);
% fprintf('临界 POF (10^-7):   %.4e\n', 1e-7);
if POF > 1e-7
    fprintf('⚠️  警告: POF 超过临界值！\n');
else
    fprintf('✓  POF 在安全范围内\n');
end
fprintf('======================================\n');

end
