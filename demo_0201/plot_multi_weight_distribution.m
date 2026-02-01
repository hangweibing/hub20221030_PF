function plot_multi_weight_distribution(xparticles, multi_weights, j, time_val, obs_a, param_obs, true_params)
%% PLOT_MULTI_WEIGHT_DISTRIBUTION 绘制各个维度的似然权重分布对比图
%
% 输入:
%   xparticles    - 当前步的粒子状态 [N x 5] (a, log_theta1, theta2, theta3, k2)
%   multi_weights - 对应的似然权重矩阵 [N x 5]
%   j             - 当前观测步编号
%   time_val      - 当前观测时间 (h)
%   obs_a         - 裂纹长度观测值 (z(j))
%   param_obs     - 参数观测值向量 [log_theta1, theta2, theta3, k2] (来自拟合)
%   true_params   - 真实参数向量 [log_theta1, theta2, theta3, k2]

    N = size(xparticles, 1);
    labels = {'Crack Length (a)', 'log-theta1', 'theta2', 'theta3', 'k2'};
    units = {'mm', '-', '-', '-', '-'};
    
    % 设置颜色和字体
    color_main = [0 0.4470 0.7410]; % 酷蓝色
    color_obs = [1 0 0];           % 红色 (观测)
    color_true = [0 0.7 0];         % 绿色 (真实)
    font_name = 'Times New Roman';
    
    figure('Name', sprintf('Observation %d (Time: %.2f h) - Multi-Weight Distribution', j, time_val), ...
           'Position', [100, 100, 1400, 800], 'Color', 'w');
    
    for d = 1:5
        subplot(2, 3, d);
        
        X = xparticles(:, d);
        W = multi_weights(:, d);
        
        % 绘制散点图 (横坐标为值，纵坐标为似然权重)
        scatter(X, W, 15, 'filled', 'MarkerFaceColor', color_main, 'MarkerFaceAlpha', 0.4);
        hold on;
        
        % 确定对应的观测值和真实值
        if d == 1
            obs_val = obs_a;
            t_val = []; % 裂纹长度真实值暂时不传，除非后续需要
        else
            if nargin >= 6 && ~isempty(param_obs)
                obs_val = param_obs(d-1);
            else
                obs_val = [];
            end
            
            if nargin >= 7 && ~isempty(true_params)
                t_val = true_params(d-1);
            else
                t_val = [];
            end
        end
        
        y_lim = get(gca, 'YLim');
        h_lgd = [];
        lgd_str = {'Particles'};
        
        % 画拟合得到的观测值 (红色虚线)
        if ~isempty(obs_val)
            h1 = plot([obs_val, obs_val], y_lim, 'r--', 'LineWidth', 1.5);
            h_lgd = [h_lgd, h1];
            lgd_str{end+1} = 'Observation (Fit)';
        end
        
        % 画真实值 (绿色实线)
        if ~isempty(t_val)
            h2 = plot([t_val, t_val], y_lim, 'Color', color_true, 'LineStyle', '-', 'LineWidth', 2);
            h_lgd = [h_lgd, h2];
            lgd_str{end+1} = 'True Value';
        end
        
        if d == 1 && ~isempty(obs_val)
             lgd_str{2} = 'Observation (Inspected)';
        end

        legend([gca().Children(end), h_lgd], lgd_str, 'Location', 'best', 'FontSize', 8);
        
        % 修饰图像
        title(sprintf('%s Likelihood Index', labels{d}), 'FontName', font_name, 'FontSize', 12, 'FontWeight', 'bold');
        xlabel(sprintf('Value [%s]', units{d}), 'FontName', font_name, 'FontSize', 10);
        ylabel('Likelihood Weight', 'FontName', font_name, 'FontSize', 10);
        grid on;
        set(gca, 'FontName', font_name, 'FontSize', 10);
    end
    
    % 在页面底部添加总标题
    sgtitle(sprintf('Multi-Weight Particle Distribution at Observation step %d', j), ...
            'FontSize', 14, 'FontWeight', 'bold', 'FontName', font_name);

    % --- 保存结果图 (新增) ---
    save_dir = 'MultiWeight_Dist_Results';
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    filename = sprintf('%s/Obs_%02d_Distribution.png', save_dir, j);
    saveas(gcf, filename);
    fprintf('  [可视化] 已保存权重分布图: %s\n', filename);
end
