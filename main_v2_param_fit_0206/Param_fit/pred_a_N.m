function [mse, a_pred] = pred_a_N(p,f,k,Const_pair_now)
% FIT_A_N 根据模型和参数，预测裂纹长度并计算 MSE
%
% 输入:
%   p       - 包含模型参数 (theta, f, k) 和代理模型 eq_fun 的结构体
%             以及几何参数 (W, B, cyclesperhour)
%   t_check - 检查时间点数组 [t1, t2, t3, ...] (小时)
%   z       - 检查得到的裂纹长度数组 [a1, a2, a3, ...] (mm)
% t_check = [4306.6400  8613.2801 12919.9201 17226.5601 21533.2001 25839.8402 30146.4802 34453.1202 38759.7603];
% z = [11.0380 11.5132 11.8109 12.4910 13.5900 15.6101 17.2276 21.5420 30.4055];
%真实参数 (log_theta1, theta2, theta3, k2): 
% [-12.5479   3.9835   0.3326   6.6750]
% 输出:
%   mse     - 预测值与检查值之间的均方误差
%   a_pred  - 预测的裂纹长度数组

    %% (1) & (2) 加载载荷谱并转换检查点为循环数
    % 加载完整载荷谱 (参考 extract_sif_cycles.m#L13-14)
    t_check=p.data_a_N(:,1);
    z=p.data_a_N(:,2);
    spectra = p.spectra; % 直接使用主程序中处理好的载荷谱
    
    % 获取循环数转换参数 
    if isfield(p, 'cyclesperhour')
        cyclesperhour = p.cyclesperhour;
    else
        cyclesperhour = 1950.70866; 
    end
    
    % 计算检查点对应的循环数 N1, N2, N3...
    N_check = round(t_check * cyclesperhour);
    
    %% 初始化迭代
    a_pred = zeros(size(z));
    a_pred(1) = z(1); % 以 (a1, N1) 作为起始点
    a_now = z(1);
    
    % 从结构体 p 获取关键参数 (与主循环保持一致)
    if isfield(p, 'step_size'), step_size = p.step_size; else, step_size = 1000; end
    if isfield(p, 'spectrum_factor'), spectrum_factor = p.spectrum_factor; else, spectrum_factor = 60; end
    if isfield(p, 'W'), W = p.W; else, W = 60; end
    if isfield(p, 'B'), B = p.B; else, B = 5; end
    if isfield(p, 'ref_load'), ref_load = p.ref_load; else, ref_load = 100; end
    
    %% (3) 迭代流程 (连续轨迹仿真)
    % N_check(1) 是起始点，仿真直到 N_check(end)
    N_now = N_check(1);
    N_total_end = N_check(end);
    
    % 当前等待记录的检查点索引
    next_check_idx = 2;
    
    while N_now < N_total_end
        % 计算步长：要么是标准 step_size，要么是到达下一个检查点的距离
        actual_step = min(step_size, N_check(next_check_idx) - N_now);
        
        % 计算当前裂纹长度下的基准 K
        K_base = sim_K_func(a_now, ref_load, W, B);
        if ~isreal(K_base)
           error('66');
        end
        % 提取载荷循环
        idx_start = 2 * N_now + 1;
        idx_end = 2 * (N_now + actual_step);
        
        if idx_end > length(spectra)
            break;
        end
        
        loads_segment = spectra(idx_start : idx_end);
        Smin = loads_segment(1:2:end);
        Smax = loads_segment(2:2:end);
        
        % 将应力谱转换为外载荷谱 (Stress to Force)
        Pmin = Smin * spectrum_factor;
        Pmax = Smax * spectrum_factor;
        
        % 根据比例关系得到 deltaK 和 Kmax 
        delta_K = (Pmax - Pmin) / ref_load * K_base;
        Kmax = Pmax / ref_load * K_base;
        
        %% (4) 计算 da 增量 (m)
        dk_m = delta_K / sqrt(1000);
        kmax_m = Kmax / sqrt(1000);
        R_ratio=1-dk_m./kmax_m;
        delta_K_th=k(2).*(1-R_ratio);
        idx_no_growth = (dk_m-delta_K_th) < 0;  % 找到 test_temp 中小于 0 的位置  
        % 计算增量
        delta_a_singgle_N = p.eq_fun(dk_m, Const_pair_now, f, k, kmax_m);
        delta_a_singgle_N(idx_no_growth) = 0;    
        delta_a_singgle_N=10.^(delta_a_singgle_N);
        %% (5) 更新裂纹长度并推进时钟
        a_now = a_now + sum(delta_a_singgle_N) * 1e3;
        if a_now>59
           a_now=59;
        end        
        if ~isreal(a_now)
           error('66');
        end
        N_now = N_now + actual_step;
        
        % 如果正好到达或越过当前的检查点，记录预测值
        if N_now >= N_check(next_check_idx)
            a_pred(next_check_idx) = a_now;
            next_check_idx = next_check_idx + 1;
            if next_check_idx > length(N_check)
                break;
            end
        end
    end
    
    %% (6) 计算 MSE
    % 只有在 z 为有效数值时计算
    mse = mean((a_pred - z).^2) / 2;
    
end
