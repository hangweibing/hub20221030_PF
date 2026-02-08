function plot_asteix_spectra()
    %% PLOT_ASTEIX_SPECTRA 绘制 AsteixSpectraData.mat 的载荷谱折线图
    %
    % 该函数会自动加载当前目录下的 AsteixSpectraData.mat 文件，
    % 并以美观的高质量图表展示其载荷循环序列。
    %
    % 使用方法:
    %   plot_asteix_spectra()

    %% 1. 数据加载与预处理
    FILENAME = 'AsteixSpectraData.mat';
    
    if ~exist(FILENAME, 'file')
        uiwait(msgbox(['未找到文件: ', FILENAME], '错误', 'error'));
        error('文件 %s 未找到，请确保该文件在 MATLAB 当前工作路径下。', FILENAME);
    end
    
    % 加载数据
    vars = load(FILENAME);
    if isfield(vars, 'spectra')
        spectra = vars.spectra;
    else
        error('变量 "spectra" 未在 %s 中找到。', FILENAME);
    end
    
    % 确保为行向量
    spectra = reshape(spectra, 1, []);
    
    %% 2. 图形参数设置 (Premium Aesthetics)
    % 定义配色方案
    color_main = [0.00, 0.45, 0.74];  % 深蓝色
    color_acc  = [0.85, 0.33, 0.10];  % 橙色 (用于高亮或缩略图)
    
    % 创建画布
    fig = figure('Name', 'Load Spectrum Visualization', ...
                 'Color', 'w', ...
                 'Units', 'normalized', ...
                 'Position', [0.1, 0.2, 0.8, 0.6]);
    
    %% 3. 绘制主图 (完整序列)
    subplot(2, 1, 1);
    h_main = plot(spectra, 'Color', [color_main 0.7], 'LineWidth', 0.5);
    hold on;
    
    % 添加标题和标签
    title('Asteix Load Spectrum - Full Sequence', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
    ylabel('Stress / Force Value', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Sequence Point Index', 'FontSize', 12);
    
    % 优化轴表现
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'Box', 'off', 'TickDir', 'out');
    set(gca, 'XLim', [1, length(spectra)]);
    
    %% 4. 绘制细节图 (展示前 200 个点以观察循环特征)
    subplot(2, 1, 2);
    detail_range = 1:min(200, length(spectra));
    plot(detail_range, spectra(detail_range), '-o', ...
         'Color', color_main, ...
         'LineWidth', 1.5, ...
         'MarkerSize', 4, ...
         'MarkerFaceColor', color_acc, ...
         'MarkerEdgeColor', 'none');
    
    title('Detail View (First ~100 Cycles)', 'FontSize', 14, 'FontWeight', 'normal');
    ylabel('Stress Value', 'FontSize', 11);
    xlabel('Sequence Point Index', 'FontSize', 11);
    
    % 添加图例说明
    legend('Load Curve', 'Peaks/Valleys', 'Location', 'best');
    
    % 优化轴表现
    grid on;
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.4, 'Box', 'off', 'TickDir', 'out');
    
    %% 5. 交互式提示 (控制台输出)
    fprintf('--------------------------------------------------\n');
    fprintf('载荷谱特征提取成功:\n');
    fprintf('  - 文件名称: %s\n', FILENAME);
    fprintf('  - 数据点总数: %d\n', length(spectra));
    fprintf('  - 峰值最大值: %.4f\n', max(spectra));
    fprintf('  - 谷值最小值: %.4f\n', min(spectra));
    fprintf('图表已生成，您可以使用放大工具观察细节。\n');
    fprintf('--------------------------------------------------\n');

end
