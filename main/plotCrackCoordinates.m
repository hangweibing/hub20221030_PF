function plotCrackCoordinates(yNewSet, zNewSet, stage_name)
%PLOTCRACKCOORDINATES 绘制裂纹坐标散点图并保存
%
% 输入参数：
%   yNewSet    - 裂纹y坐标数组
%   zNewSet    - 裂纹z坐标数组
%   stage_name - 阶段名称（用于标题和文件名）
%
% 功能：
%   绘制裂纹节点的散点图，显示物理边界，并保存为SVG矢量图

    % 创建图形（不显示，后台保存）
    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600]);
    
    % 绘制裂纹节点
    plot(real(zNewSet), real(yNewSet), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5);
    hold on;
    
    % 标注末端点
    plot(real(zNewSet(end)), real(yNewSet(end)), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
    
    % 绘制物理边界
    z_min = min(real(zNewSet)) - 2;
    z_max = max(real(zNewSet)) + 2;
    z_boundary = z_min:0.1:z_max;
    y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
    y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
    plot(z_boundary, y_upper, 'k--', 'LineWidth', 1);
    plot(z_boundary, y_lower, 'k--', 'LineWidth', 1);
    
    % 图形设置
    xlim([z_min, z_max]);
    ylim([6, 15]);
    xlabel('z坐标 (mm)', 'FontSize', 12);
    ylabel('y坐标 (mm)', 'FontSize', 12);
    title(sprintf('裂纹坐标可视化 - %s\n节点数: %d', stage_name, length(yNewSet)), 'FontSize', 14);
    grid on;
    legend('裂纹节点', '末端点', '上边界', '下边界', 'Location', 'best');
    

    
    % 保存图形
    result_folder = 'addConstraintNewSatgeFunc_debug结果';
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    filename = sprintf('%s/裂纹坐标_%s_%s.svg', result_folder, strrep(stage_name, ' ', '_'), timestamp);
    saveas(fig, filename);
    close(fig);
    
    fprintf('  [plotCrackCoordinates] 已保存图形: %s\n', filename);
end
