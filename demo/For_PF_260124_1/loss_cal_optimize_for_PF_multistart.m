function [loss_out, delta_kth_out, kc_out, theta_out, opti_mode_used, constraint_wrong, pass_index] = loss_cal_optimize_for_PF_multistart(p, pop_now, N_starts)
    % Multi-start version of loss_cal_optimize_for_PF using parfor
    % p: parameter structure (parameter_gp)
    % pop_now: cell array containing population data
    % N_starts: Number of starting points to use
    
    PARA0_default = pop_now{1,6};
    pop_all = pop_now{1,1};
    
    % Initialize starting points array
    % Based on PARA0 length (should be 5 based on user description)
    len_para = length(PARA0_default);
    all_PARA0 = zeros(N_starts, len_para);
    all_PARA0(1, :) = PARA0_default;
    
    if N_starts > 1
        total_pop = size(pop_all, 1);
        % Randomly select rows
        % Note: randperm ensures no duplicate rows are picked from pop_all
        rand_rows = randperm(total_pop, min(N_starts - 1, total_pop));
        for i = 1:length(rand_rows)
            row_idx = rand_rows(i);
            % Column mapping: k1=col 2, k2=col 3, f1=col 4, f2=col 5, f3=col 6
            % We take 5 parameters from these columns
            candidate = [pop_all(row_idx, 2), pop_all(row_idx, 3), ...
                         pop_all(row_idx, 4), pop_all(row_idx, 5), ...
                         pop_all(row_idx, 6)];
            
            % Ensure candidate matches PARA0 length if it differs for some reason
            if length(candidate) == len_para
                all_PARA0(i+1, :) = candidate;
            else
                % If mismatch, just repeat default or handle differently
                all_PARA0(i+1, :) = PARA0_default;
            end
        end
    end
    
    % Storage for results
    all_results = cell(N_starts, 1);
    
    parfor i = 1:N_starts
    % parfor i = 1:N_starts
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
            if pass == 1
                res.total_loss = l(1); 
                fprintf('  第 %d 次尝试 (Start): MSE Loss = %.4e, Parameters: [dk=%.4f, kc=%.4f, theta=[%.4f, %.4f, %.4f]]\n', ...
                        i, l(1), dk, kcc, th(1), th(2), th(3));
            else
                res.total_loss = inf;
                fprintf('  第 %d 次尝试 (Start): 拟合失败 (pass=%d, constraint_wrong=%d)\n', i, pass, cw);
            end
        catch ME
            res = struct('total_loss', inf, 'pass_index', 0, 'loss', [inf, inf, inf], ...
                         'delta_kth', 0, 'kc', 0, 'theta', zeros(3,1), 'opti_mode_used', 0, 'constraint_wrong', 0);
            fprintf('  第 %d 次尝试 (Start): 异常退出 (%s)\n', i, ME.message);
        end
        all_results{i} = res;
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
