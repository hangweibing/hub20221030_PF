function [loss_out, delta_kth_out, kc_out, theta_out, opti_mode_used, constraint_wrong, pass_index] = loss_cal_optimize_for_PF_multistart(p, pop_now, N_starts, historical_PARA0)
% Multi-start version of loss_cal_optimize_for_PF using parfor
% p: parameter structure (parameter_gp)
% pop_now: cell array containing population data
% N_starts: Number of starting points to use
% historical_PARA0: (Optional) Matrix of parameters from all previous fitting steps [num_prev x len_para]

PARA0_default = pop_now{1,6};
pop_all = pop_now{1,1};

% Initialize starting points array
len_para = length(PARA0_default);
if isfield(p, 'fix_k2') && p.fix_k2
    PARA0_default(2) = p.k2_fixed_val;
end
if isfield(p, 'fix_theta3') && p.fix_theta3
    PARA0_default(4) = p.theta3_fixed_val;
end
all_PARA0 = zeros(N_starts, len_para);
all_PARA0(1, :) = PARA0_default;

start_idx = 2;
% 如果提供了历史拟合结果，将其作为后续起始点
if nargin >= 4 && ~isempty(historical_PARA0)
    [num_hist, col_hist] = size(historical_PARA0);
    if col_hist == len_para
        % 最多取 N_starts-1 个历史点，防止超出总起始点数
        num_to_use = min(num_hist, N_starts - 1);
        all_PARA0(2 : 1 + num_to_use, :) = historical_PARA0(end - num_to_use + 1 : end, :);
        start_idx = 2 + num_to_use;
        fprintf('  [Multistart] 将 %d 组历史拟合结果作为起始点。\n', num_to_use);
    end
end

if N_starts >= start_idx
    total_pop = size(pop_all, 1);
    num_to_sample = N_starts - (start_idx - 1);
    % 从 population 中随机采样剩余的起始点
    rand_indices = randperm(total_pop, num_to_sample);

    for i = 1:num_to_sample
        row_idx = rand_indices(i);
        % Column mapping: k1=col 2, k2=col 3, f1=col 4, f3=col 5, f2=col 6
        k1_val = PARA0_default(1); 
        k2_val = pop_all(row_idx, 3);
        if isfield(p, 'fix_k2') && p.fix_k2
            k2_val = p.k2_fixed_val;
        end
        f1_val = pop_all(row_idx, 4);
        f2_val = pop_all(row_idx, 5); 
        f3_val = pop_all(row_idx, 6);
        if isfield(p, 'fix_theta3') && p.fix_theta3
            f2_val = p.theta3_fixed_val;
        end
        candidate = [k1_val, k2_val, f1_val, f2_val, f3_val];
        all_PARA0(start_idx + i - 1, :) = candidate;
    end
end

% Storage for results
all_results = cell(N_starts, 1);

% parfor i = 1:N_starts
parfor i = 1:N_starts
    p_local = p;
    % Set current starting point
    p_local.PARA0 = all_PARA0(i, :);
    p_local.parameter_K = p_local.PARA0(1:2);
    p_local.theta = p_local.PARA0(3:end);

    % Call original optimization function
    try
        [l, dk, kcc, th, mode, cw, pass] = loss_cal_optimize_for_PF(p_local);

        res = struct();
        res.loss = l;
        res.delta_kth = dk;
        res.kc = kcc;
        res.theta = th;
        res.opti_mode_used = mode;
        res.constraint_wrong = cw;
        res.pass_index = pass;

        % 只把参数拟合后的MSE作为loss进行选择
        res.total_loss = l(1);
        fprintf('  第 %d 次尝试 (Start): Loss=[MSE:%.4e, P1:%.4e, P2:%.4e], Parameters: [dk=%.4f, kc=%.4f, theta=[%.4f, %.4f, %.4f]]\n', ...
            i, l(1), l(2), l(3), dk, kcc, th(1), th(2), th(3));
    catch ME
        res = struct('total_loss', inf, 'pass_index', 0, 'loss', [inf, inf, inf], ...
            'delta_kth', 0, 'kc', 0, 'theta', zeros(3,1), 'opti_mode_used', 0, 'constraint_wrong', 0);
        fprintf('  第 %d 次尝试 (Start): 异常退出 (%s)\n', i, ME.message);
    end
    all_results{i} = res;
end

% --- 结果处理与方差计算 ---
fprintf('\n[多起始点优化结果汇总]\n');
valid_params = [];
valid_mse = [];

for i = 1:N_starts
    r = all_results{i};
    if r.pass_index == 1
        % 记录 4 个主要优化参数: theta (3个) + kc (1个)
        % 注意: 如果 fix_k2=false，则 kc 也是变量
        params_row = [r.theta', r.kc];
        valid_params = [valid_params; params_row];
        valid_mse = [valid_mse; r.total_loss];
    end
end

% 提取 MSE < 1e-3 的组合
threshold = 1e-3;
idx_filter = valid_mse < threshold;
filtered_params = valid_params(idx_filter, :);
filtered_mse = valid_mse(idx_filter);

if ~isempty(filtered_params)
    fprintf('找到 %d 组 MSE < %.0e 的参数组合:\n', size(filtered_params, 1), threshold);
    % 计算方差
    p_vars = var(filtered_params, 1); % 使用 N 归一化
    fprintf('  参数方差 (f1, f3, f2, kc): [%.4e, %.4e, %.4e, %.4e]\n', ...
        p_vars(1), p_vars(2), p_vars(3), p_vars(4));

    % 如果有多组，打印第一组和最后一组示例
    if size(filtered_params, 1) > 1
        fprintf('  示例组合 1: MSE=%.4e, Params=[%.4f, %.4f, %.4f, %.4f]\n', ...
            filtered_mse(1), filtered_params(1,1), filtered_params(1,2), filtered_params(1,3), filtered_params(1,4));
    end
else
    fprintf('未找到 MSE < %.0e 的参数组合。\n', threshold);
end

% Find the best result based on minimal total loss
total_losses = cellfun(@(x) x.total_loss, all_results);
[min_loss, best_idx] = min(total_losses);

if min_loss < inf
    best_res = all_results{best_idx};
else
    % If none succeeded, return results from the first start (default start)
    best_res = all_results{1};
end

% Assign outputs
loss_out = best_res.loss(1); % 只返回MSE
delta_kth_out = best_res.delta_kth;
kc_out = best_res.kc;
theta_out = best_res.theta;
opti_mode_used = best_res.opti_mode_used;
constraint_wrong = best_res.constraint_wrong;
pass_index = best_res.pass_index;

if min_loss < inf
    fprintf('  [Multistart完成] 最佳尝试: %d, 最小MSE: %.4e\n', best_idx, min_loss);
else
    fprintf('  [Multistart失败] 所有尝试均未成功拟合\n');
end
end
