%% ===================================================================
%% 程序名称：gen_synthetic_data_hub.m
%% 功能描述：参考 hub20260127.m 的循环逻辑生成合成观测数据
%%           用于主程序的检查数据生成，支持手动设置参数行数或随机选取
%% ===================================================================

% 清理工作空间
clear; close all; clc;

% 设置随机数种子 (确保可重复性)
SIM_SEED = 2026;
rng(SIM_SEED);

%% 1. 环境配置与参数加载
fprintf('--- 环境配置与参数加载 ---\n');

% 加载载荷谱 (参考 hub20260127.m#L289-295)
if exist('AsteixSpectraData.mat', 'file')
    load('AsteixSpectraData.mat', 'spectra');
else
    error('未找到 AsteixSpectraData.mat，请确保在正确的目录下运行。');
end
% 重复载荷谱以确保足够长
spectra = repmat(reshape(spectra, 1, []), 1, 100);
spectra = spectra(2:end);

% 加载 POD 模型数据 (包含 Uinput_integrated, Uinput_splitted 等核心模型参数)
if exist('pod_models.mat', 'file')
    load('pod_models.mat');
else
    error('未找到 pod_models.mat，请确保在正确的目录下运行。');
end

% 加载真实参数库 (AM-TC4-GRO.xlsx)
if exist('AM-TC4-GRO.xlsx', 'file')
    T = readtable('AM-TC4-GRO.xlsx');
    params_pool = table2array(T);
    num_data = size(params_pool, 1);
else
    error('未找到 AM-TC4-GRO.xlsx，请检查路径。');
end

% 阶段映射与测试误差 (参考 hub20260127.m#L352-358)
stage_map = [1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 8];
testErrSet = [0.019803420755871 0.058377623733943 0.013343080193008 ...
    0.009221345729625 0.037283991994545 0.011408674471790 ...
    0.082711915490345 0.032375406062069 0.041104885753232 ...
    0.057544490601307];

% =========================================================================
% 参数选取设置
% =========================================================================
% 您可以在此处手动设置选取的行数，或者保持随机选取
manual_idx = 17; % 修改此处手动指定行数 (17)
if isempty(manual_idx)
    idx = randi(num_data);
else
    idx = manual_idx;
end

true_params = params_pool(idx, :);
log_theta1_ = true_params(1);
theta2 = true_params(2);
theta3 = true_params(3);
k2 = true_params(4);

fprintf('选取参数行数: %d / %d\n', idx, num_data);
fprintf('真实参数组合: log_theta1=%.4f, theta2=%.4f, theta3=%.4f, k2=%.4f\n', ...
    log_theta1_, theta2, theta3, k2);

%% 2. 仿真参数与初始状态初始化
fprintf('--- 仿真初始化 ---\n');

% 裂纹几何参数 (参考 hub20260127.m#L207-257)
n_nodes = 21;
centers = [13, 30];
thetas = linspace(3/2*pi, 2*pi, n_nodes);

a_init = 2.0; % 初始裂纹半径 (mm)
c_init = a_init;
yRegSet = centers(1) + c_init * sin(thetas);
zRegSet = centers(2) + a_init * cos(thetas);

% 其他仿真常数
step = 1000;              % 时间步长 (cycles)
cyclesperhour = 1950.70866; % 每小时循环数
SPLITTED = 0;             % 初始分裂状态

% 分裂模型数据映射
Uinput_splitted_1 = Uinput_splitted{1};
averInput_splitted_1 = averInput_splitted{1};
Uinput_splitted_2 = Uinput_splitted{2};
averInput_splitted_2 = averInput_splitted{2};

% --- 轨迹记录设置 ---
DEBUG_MODE = true;              % 是否记录轨迹
record_interval = 30;           % 每隔多少个循环记录一次
particles_coordinates_history = []; % 用于记录 [y, z] 坐标

%% 3. 载荷谱平均增量计算 (参考 hub20260127.m#L300-345)
num_steps = floor((length(spectra))/2/step);
aver_delta_sigma_set = zeros(1, num_steps);
aver_R_set = zeros(1, num_steps);
aver_Smax_set = zeros(1, num_steps);

k_ptr = 1;
for i = 1:num_steps
    idx_seg = 2*(k_ptr + (1:step) - 1);
    Smax_vec = spectra(idx_seg);
    Smin_vec = spectra(idx_seg - 1);

    aver_delta_sigma_set(i) = mean(Smax_vec - Smin_vec);
    aver_R_set(i) = mean(Smin_vec) / mean(Smax_vec);
    aver_Smax_set(i) = mean(Smax_vec);
    k_ptr = k_ptr + step;
end

%% 4. 裂纹扩展主循环
fprintf('--- 开始裂纹扩展仿真 ---\n');

t_history = [0];
z_history = [zRegSet(end)]; % 记录上表面 z 坐标
m = 1;

while zRegSet(end) < 50
    if m > num_steps
        fprintf('警告: 载荷谱已耗尽，z 坐标未达到 50。\n');
        break;
    end

    % 获取当前裂纹尺寸用于阶段判断
    a_up = zRegSet(end) - 30;
    a_down = zRegSet(1) - 30;

    % 获取阶段索引 (参考 hub20260127.m#L446)
    m_index = getModelIndexFunc(a_up, a_down);
    if m_index == 0
        fprintf('错误: getModelIndexFunc 返回 0 (无效状态)，仿真终止。\n');
        break;
    end

    % 模型选择 (参考 hub20260127.m#L457-470)
    if SPLITTED && (m_index==3 || m_index==5)
        if m_index==3
            curUinput = Uinput_splitted_1;
            curAverInput = averInput_splitted_1;
        elseif m_index==5
            curUinput = Uinput_splitted_2;
            curAverInput = averInput_splitted_2;
        end
        m_name = sprintf('nn_stage%ds', stage_map(m_index));
    else
        curUinput = Uinput_integrated{stage_map(m_index)};
        curAverInput = averInput_integrated{stage_map(m_index)};
        m_name = sprintf('nn_stage%d', stage_map(m_index));
    end

    % 执行一步扩展 (使用 a2aNew 包装函数)
    % 注: 此处作为 "True" 数据生成，不显式添加额外噪声（a2aNew 内部可能包含固定的小随机性或 POD 误差）
    [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, ~] = ...
        a2aNew(yRegSet, zRegSet, aver_delta_sigma_set(m), aver_R_set(m), aver_Smax_set(m), m_name, ...
        curUinput, curAverInput, log_theta1_, theta2, theta3, k2, step, testErrSet, 1);

    % 记录历史
    m = m + 1;
    current_hour = (m-1) * step / cyclesperhour;
    t_history(end+1) = current_hour;
    z_history(end+1) = zRegSet(end);

    % 记录坐标轨迹
    if DEBUG_MODE && mod(m-1, record_interval) == 0
        particles_coordinates_history = [particles_coordinates_history; [yRegSet(:)', zRegSet(:)']];
    end

    % 打印进度
    if mod(m, 100) == 0
        fprintf('  当前时间: %.2f h, 裂纹 z 坐标: %.4f mm\n', current_hour, zRegSet(end));
    end
end

total_life = t_history(end);
fprintf('仿真结束。总寿命: %.2f h, 最终 z 坐标: %.4f mm\n', total_life, zRegSet(end));

%% 5. 等间隔记录 10 组检查数据
num_checks = 10;
t_check = linspace(100, total_life, num_checks);

% 插值得到对应时间的 z 坐标 (合成观测值)
z_inspected = interp1(t_history, z_history, t_check, 'linear');

% 打印结果以供直接复制到 hub20260127.m
fprintf('\n===========================================\n');
fprintf('   生成的检查数据 (用于 hub20260127.m)\n');
fprintf('===========================================\n');
fprintf('t_check = [%s];\n', num2str(t_check, '%.4e '));
fprintf('z       = [%s];\n', num2str(z_inspected, '%.4e '));
fprintf('-------------------------------------------\n');
fprintf('真实参数组合: [%s]\n', num2str(true_params, '%.4e '));
fprintf('===========================================\n');

%% 6. 绘图显示
figure('Name', 'Synthetic Data Generation (HUB)', 'Color', 'w', 'Position', [200 200 800 500]);
plot(t_history, z_history, 'b-', 'LineWidth', 2, 'DisplayName', 'True Growth Path');
hold on;
plot(t_check, z_inspected, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Inspections');
xlabel('Flight Hours (h)');
ylabel('Crack Z-coordinate (mm)');
title(['Synthetic Data (Param Row: ', num2str(idx), ')']);
grid on;
legend('Location', 'best');

% 保存结果到 mat file
save('synthetic_data_hub.mat', 't_check', 'z_inspected', 'true_params', 'idx', 't_history', 'z_history');
if DEBUG_MODE
    save('synthetic_data_hub.mat', 'particles_coordinates_history', '-append');
end
fprintf('数据已保存至 synthetic_data_hub.mat\n');

%% 7. 绘制裂纹轨迹图 (参考 debug_plot_particle_trajectories.m)
if DEBUG_MODE && ~isempty(particles_coordinates_history)
    particle_data = particles_coordinates_history;
    [num_time_steps, num_coords] = size(particle_data);
    num_nodes = num_coords / 2;

    % 创建图形
    fig_traj = figure('Name', 'Crack Trajectory', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
    colors = lines(num_time_steps);
    hold on;

    % --- 绘制物理边界 ---
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 区域1: z从30~35，y=13和y=7两条平行线
    z_region1 = 30:0.1:35;
    y_upper_region1 = 13 * ones(size(z_region1));
    y_lower_region1 = 7 * ones(size(z_region1));

    % 区域2: z从35到37.82842712，通过两个圆弧连接
    z_region2 = 35:0.01:37.82842712;
    % 上边界圆弧
    y_upper_region2 = upCenter(2) - sqrt(radius^2 - (z_region2 - upCenter(1)).^2);
    % 下边界圆弧
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
    plot(z_boundary, y_upper, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Upper Boundary');
    plot(z_boundary, y_lower, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Lower Boundary');

    % --- 绘制裂纹轨迹 ---
    for i_step = 1:num_time_steps
        y_c = particle_data(i_step, 1:num_nodes);
        z_c = particle_data(i_step, num_nodes+1:end);

        plot(real(z_c), real(y_c), '-', 'Color', colors(i_step, :), 'LineWidth', 1.5);
        % 最后一个步额外突出显示
        if i_step == num_time_steps
            plot(real(z_c), real(y_c), 'ro-', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Current Crack');
        end
    end

    % 图形设置
    xlim([30, 52]);
    ylim([5, 15]);
    axis equal;
    xlabel('z Coordinate (mm)');
    ylabel('y Coordinate (mm)');
    title(sprintf('Crack Growth Trajectory (Interval: %d steps)', record_interval));
    grid on;
    % 只显示边界和最后一帧的图例，避免太拥挤
    legend('show', 'Location', 'best');
end
