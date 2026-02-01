function visualize_parameter_distribution(particle_params, training_data_file)
%% VISUALIZE_PARAMETER_DISTRIBUTION 可视化采样参数与实验数据的对比
%
% 输入:
%   particle_params    - 采样得到的粒子参数矩阵 [N x 4]
%   training_data_file - 原始实验数据文件名 (默认 'AM-TC4-GRO.xlsx')

    if nargin < 2
        training_data_file = 'AM-TC4-GRO.xlsx';
    end

    % 检查文件是否存在
    if ~exist(training_data_file, 'file')
        warning('未找到实验数据文件 %s，无法进行对比。', training_data_file);
        return;
    end

    % 获取训练数据用于可视化对比
    T = readtable(training_data_file);
    X_train = table2array(T);
    param_names = T.Properties.VariableNames;
    d = size(X_train, 2);

    % 颜色配置 (Premium 风格)
    color_sampled = [0 0.4470 0.7410]; % 酷蓝色 (Sampled)
    color_train = [0.8500 0.3250 0.0980];   % 砖红色 (Training)
    font_name = 'Times New Roman';

    %% 1. 单参数边缘分布对比 (直方图)
    figure('Name', '采样参数分布对比 - Histograms', 'Position', [100, 100, 1200, 700], 'Color', 'w');
    for i = 1:min(d, 4)
        subplot(2,2,i)
        
        % 绘制采样数据的分布
        histogram(particle_params(:,i), 40, 'Normalization','pdf', ...
            'FaceColor', color_sampled, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
        hold on
        
        % 绘制训练数据的分布
        histogram(X_train(:,i), 15, 'Normalization','pdf', ...
            'FaceColor', color_train, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

        % 设置标题和标签
        clean_name = strrep(param_names{i}, '_', '-');
        title(sprintf('%s Distribution Comparison', clean_name), 'Interpreter', 'none', ...
            'FontSize', 12, 'FontName', font_name, 'FontWeight', 'bold')
        xlabel(clean_name, 'Interpreter', 'none', 'FontSize', 10, 'FontName', font_name)
        ylabel('Probability Density', 'FontSize', 10, 'FontName', font_name)
        legend('Sampled (KDE)', 'Experimental Data', 'Location', 'best', ...
            'FontSize', 9, 'FontName', font_name)
        
        grid on
        set(gca, 'FontName', font_name, 'FontSize', 10, 'Box', 'off')
    end

    %% 2. 三参数联合分布可视化 (3D 散点图)
    figure('Name', '三参数联合分布对比 - 3D Scatter', 'Position', [150, 150, 1300, 850], 'Color', 'w');

    % 计算需要显示的三参数组合 (4选3)
    param_triplets = [
        1,2,3;  % log_theta1, theta2, theta3
        1,2,4;  % log_theta1, theta2, k2
        1,3,4;  % log_theta1, theta3, k2
        2,3,4   % theta2, theta3, k2
    ];

    for k = 1:size(param_triplets, 1)
        i = param_triplets(k,1);
        j = param_triplets(k,2);
        l = param_triplets(k,3);

        subplot(2,2,k);

        % 绘制实验数据的三维散点图 (红色)
        scatter3(X_train(:,i), X_train(:,j), X_train(:,l), 35, 'filled', ...
            'MarkerFaceColor', color_train, 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
        hold on;

        % 绘制采样数据的三维散点图 (蓝色)
        scatter3(particle_params(:,i), particle_params(:,j), particle_params(:,l), 20, 'filled', ...
            'MarkerFaceColor', color_sampled, 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');

        % 清理变量名用于显示
        clean_name_i = strrep(param_names{i}, '_', '-');
        clean_name_j = strrep(param_names{j}, '_', '-');
        clean_name_l = strrep(param_names{l}, '_', '-');

        xlabel(clean_name_i, 'Interpreter', 'none', 'FontSize', 10, 'FontName', font_name);
        ylabel(clean_name_j, 'Interpreter', 'none', 'FontSize', 10, 'FontName', font_name);
        zlabel(clean_name_l, 'Interpreter', 'none', 'FontSize', 10, 'FontName', font_name);

        title(sprintf('%s - %s - %s', clean_name_i, clean_name_j, clean_name_l), ...
            'Interpreter', 'none', 'FontSize', 11, 'FontName', font_name, 'FontWeight', 'bold');

        legend('Exp Data', 'KDE Sampled', 'Location', 'best', 'FontSize', 9, 'FontName', font_name);

        grid on;
        set(gca, 'FontName', font_name, 'FontSize', 10);
        view(45, 25); % 调整视角
    end

    fprintf('--- 参数分布可视化完成！ ---\n');
end
