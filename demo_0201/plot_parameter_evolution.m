function plot_parameter_evolution(Xpf_params, P_upper, P_lower, time_axis, j, true_params)
%% PLOT_PARAMETER_EVOLUTION 绘制 4 个材料参数随时间演化的曲线图（包含 90% 置信区间）
%
% 输入:
%   Xpf_params - 参数均值历史 [4 x m] (log_theta1, theta2, theta3, k2)
%   P_upper    - 参数上边界历史 [4 x m] (95th percentile)
%   P_lower    - 参数下边界历史 [4 x m] (5th percentile)
%   time_axis  - 时间轴 [1 x m] (h)
%   j          - 当前观测步编号
%   true_params- 真实参数值 [1 x 4]

    labels = {'log-theta1', 'theta2', 'theta3', 'k2'};
    font_name = 'Times New Roman';
    color_mean = [0 0.4470 0.7410]; % 均值线颜色 (蓝色)
    color_shade = [0.1 0.6 0.9];    % 置信区间阴影颜色 (浅蓝)
    color_true = [0 0.7 0];         % 真实值颜色 (绿色)

    fig = figure('Name', sprintf('Observation %d - Parameter Evolution', j), ...
           'Position', [150, 150, 1200, 800], 'Color', 'w');

    for d = 1:4
        subplot(2, 2, d);
        hold on;
        
        % 1. 绘制 90% 置信区间 (浅色背景)
        fill_x = [time_axis, fliplr(time_axis)];
        fill_y = [P_upper(d, :), fliplr(P_lower(d, :))];
        fill(fill_x, fill_y, color_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
        % 2. 绘制均值线 (实线)
        plot(time_axis, Xpf_params(d, :), 'Color', color_mean, 'LineWidth', 2);
        
        % 3. 绘制真实值 (横向虚线)
        if nargin >= 6 && ~isempty(true_params)
            plot([time_axis(1), time_axis(end)], [true_params(d), true_params(d)], ...
                'Color', color_true, 'LineStyle', '--', 'LineWidth', 1.5);
        end
        
        % 修饰图像
        title(labels{d}, 'FontName', font_name, 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Time (Flight hours)', 'FontName', font_name, 'FontSize', 11);
        ylabel('Value', 'FontName', font_name, 'FontSize', 11);
        grid on;
        set(gca, 'FontName', font_name, 'FontSize', 11);
        
        if d == 1
            if nargin >= 6
                legend('90% CI', 'Prediction Mean', 'True Value', 'Location', 'best');
            else
                legend('90% CI', 'Prediction Mean', 'Location', 'best');
            end
        end
    end

    % 添加总标题
    sgtitle(sprintf('Evolution of Material Parameters (Obs Step: %d)', j), ...
            'FontSize', 16, 'FontWeight', 'bold', 'FontName', font_name);

    % --- 保存图像 ---
    save_dir = 'MultiWeight_Dist_Results';
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    filename = sprintf('%s/Obs_%02d_Param_Evolution.png', save_dir, j);
    saveas(gcf, filename);
    fprintf('  [可视化] 已保存参数演化图: %s\n', filename);
    
    % 可视化后建议关闭，防止内存占用过多
    % close(fig); 
end
