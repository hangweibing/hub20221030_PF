% debug绘图脚本：绘制指定粒子在不同时间步下的裂纹轨迹
% 从debug_particles_coordinates.mat读取数据，根据输入的粒子编号绘制轨迹

% 加载数据
load('debug_particles_coordinates.mat', 'particles_coordinates');

% 获取用户输入的粒子编号和间隔数
particle_idx = 1023;
interval = 2;

% 检查粒子编号是否有效
if particle_idx < 1 || particle_idx > length(particles_coordinates)
    error('粒子编号超出范围 (1-%d)', length(particles_coordinates));
end

% 获取指定粒子的坐标历史
particle_data = particles_coordinates{particle_idx};
[num_time_steps, num_coords] = size(particle_data);

if num_time_steps == 0
    error('粒子%d没有坐标数据', particle_idx);
end

% 检查坐标数量是否正确（应为42个：21个y + 21个z）
if mod(num_coords, 2) ~= 0 || num_coords == 0
    error('粒子%d的坐标数据格式错误', particle_idx);
end

num_nodes = num_coords / 2;  % 节点数量（21个）

% 创建图形
fig = figure('Visible', 'on', 'Position', [100, 100, 1200, 600]);

% 定义颜色映射
colors = lines(num_time_steps);  % 使用MATLAB的lines颜色方案

% 绘制物理边界（与plotCrackCoordinates一致）
hold on;

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

% 区域3: z从37.82842712到52，y=11和y=9两条直线
z_region3 = 37.82842712:0.1:52;
y_upper_region3 = 11 * ones(size(z_region3));
y_lower_region3 = 9 * ones(size(z_region3));

% 合并所有区域
z_boundary = [z_region1, z_region2, z_region3];
y_upper = [y_upper_region1, y_upper_region2, y_upper_region3];
y_lower = [y_lower_region1, y_lower_region2, y_lower_region3];

% 绘制边界
plot(z_boundary, y_upper, 'k-', 'LineWidth', 1.5);
plot(z_boundary, y_lower, 'k-', 'LineWidth', 1.5);

% 绘制间隔的时间步轨迹
legend_entries = {};
plot_indices = 1:interval:num_time_steps;  % 要绘制的时间步索引
num_plots = length(plot_indices);

for i = 1:num_plots
    m = plot_indices(i);  % 实际的时间步索引

    % 获取第m个时间步的坐标数据
    % 每行前num_nodes个元素是y坐标，后num_nodes个元素是z坐标
    y_coords = particle_data(m, 1:num_nodes);
    z_coords = particle_data(m, num_nodes+1:end);

    % 绘制轨迹
    plot(real(z_coords), real(y_coords), 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
         'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));

    % 添加图例条目
    legend_entries{end+1} = sprintf('时间步 %d', m);
end

% 图形设置
xlim([30, 52]);
ylim([6, 14]);
axis equal;
xlabel('z坐标 (mm)', 'FontSize', 12);
ylabel('y坐标 (mm)', 'FontSize', 12);
title(sprintf('粒子%d的裂纹轨迹 - 总%d个时间步, 绘制%d个轨迹 (间隔%d), %d个节点', ...
    particle_idx, num_time_steps, num_plots, interval, num_nodes), 'FontSize', 14);
grid on;

% 添加图例
legend_entries = [{'上边界', '下边界'}, legend_entries];
legend(legend_entries, 'Location', 'best', 'NumColumns', 2);

fprintf('已绘制粒子%d的%d个轨迹 (间隔%d, 总%d个时间步)\n', particle_idx, num_plots, interval, num_time_steps);