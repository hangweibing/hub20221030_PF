function plotCrackCoordinates(yNewSet, zNewSet, particle_idx, time_step, flight_hours)
%PLOTCRACKCOORDINATES 绘制裂纹坐标散点图并保存
%
% 输入参数：
%   yNewSet     - 裂纹y坐标数组
%   zNewSet     - 裂纹z坐标数组
%   particle_idx - 粒子编号（可选，用于标题和文件名）
%   time_step   - 时间步（可选，用于标题和文件名）
%   flight_hours - 飞行小时数（可选，用于标题和文件名）
%
% 功能：
%   绘制裂纹节点的散点图，显示物理边界，并保存为SVG矢量图

    % 创建图形（不显示，后台保存）
    % 调整窗口大小以匹配坐标轴比例
    % x轴范围: 20单位, y轴范围: 6单位, 比例为20:6=10:3≈3.33:1
    % 设置窗口宽度1200像素，则高度应为1200/(20/6)≈360像素
    fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 360]);
    
    % 绘制裂纹节点
    plot(real(zNewSet), real(yNewSet), 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5);
    hold on;
    
    % 标注末端点
    plot(real(zNewSet(end)), real(yNewSet(end)), 'r^', 'MarkerSize', 10, 'LineWidth', 2);
    
    % 绘制物理边界
    % 圆心坐标
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 区域1: z从30~35，y=13和y=7两条平行线
    z_region1 = 30:0.1:35;
    y_upper_region1 = 13 * ones(size(z_region1));
    y_lower_region1 = 7 * ones(size(z_region1));

    % 区域2: z从35到37.82842712，通过两个圆弧连接
    z_region2 = 35:0.1:37.82842712;
    % 上边界圆弧（圆心[37.82842712, 14]，半径3）
    y_upper_region2 = upCenter(2) - sqrt(radius^2 - (z_region2 - upCenter(1)).^2);
    % 下边界圆弧（圆心[37.82842712, 6]，半径3）
    y_lower_region2 = downCenter(2) + sqrt(radius^2 - (z_region2 - downCenter(1)).^2);

    % 区域3: z从37.82842712到50，y=11和y=9两条直线
    z_region3 = 37.82842712:0.1:50;
    y_upper_region3 = 11 * ones(size(z_region3));
    y_lower_region3 = 9 * ones(size(z_region3));

    % 合并所有区域
    z_boundary = [z_region1, z_region2, z_region3];
    y_upper = [y_upper_region1, y_upper_region2, y_upper_region3];
    y_lower = [y_lower_region1, y_lower_region2, y_lower_region3];

    % 绘制边界
    plot(z_boundary, y_upper, 'k-', 'LineWidth', 1);
    plot(z_boundary, y_lower, 'k-', 'LineWidth', 1);
    
    % 图形设置
    xlim([30, 50]);
    ylim([7, 13]);
    axis equal;  % 使横纵坐标单位长度相等
    xlabel('z坐标 (mm)', 'FontSize', 12);
    ylabel('y坐标 (mm)', 'FontSize', 12);
    
    % 构建标题（根据可用参数）
    if nargin >= 5 && ~isempty(flight_hours)
        title_str = sprintf('裂纹坐标可视化 - 粒子%d, 时间步%d, %.2f小时\n节点数: %d', ...
            particle_idx, time_step, flight_hours, length(yNewSet));
    elseif nargin >= 3 && ~isempty(particle_idx)
        title_str = sprintf('裂纹坐标可视化 - 粒子%d, 时间步%d\n节点数: %d', ...
            particle_idx, time_step, length(yNewSet));
    else
        title_str = sprintf('裂纹坐标可视化\n节点数: %d', length(yNewSet));
    end
    title(title_str, 'FontSize', 14);
    grid on;
    legend('裂纹节点', '末端点', '上边界', '下边界', 'Location', 'best');
    

    
    % 保存图形
    result_folder = 'addConstraintNewSatgeFunc_debug结果';
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    if nargin >= 3 && ~isempty(particle_idx)
        filename = sprintf('%s/裂纹坐标_Particle%d_Step%d_%s.svg', ...
            result_folder, particle_idx, time_step, timestamp);
    else
        filename = sprintf('%s/裂纹坐标_%s.svg', result_folder, timestamp);
    end
    saveas(fig, filename);
    close(fig);
    
    fprintf('  [plotCrackCoordinates] 已保存图形: %s\n', filename);
end
