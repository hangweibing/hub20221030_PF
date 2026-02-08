function [mse, a_pred] = pred_a_N(p, f, k, Const_pair_now)
% PRED_A_N 使用完整几何信息和神经网络模型预测裂纹长度并计算 MSE
%
% 输入:
%   p       - 包含模型参数和配置的结构体
%   f       - theta 参数向量 [theta1, theta2, theta3]
%   k       - K 参数向量 [delta_kth, kc]
%   Const_pair_now - 常数对
%
% 输出:
%   mse     - 预测值与观测值之间的均方误差
%   a_pred  - 预测的裂纹长度数组

%% (1) 提取配置参数
t_check = p.data_a_N(:, 1);
z = p.data_a_N(:, 2);
spectra = p.spectra;

% 获取循环数转换参数
if isfield(p, 'cyclesperhour')
    cyclesperhour = p.cyclesperhour;
else
    cyclesperhour = 1950.70866;
end

% 计算检查点对应的循环数
N_check = round(t_check * cyclesperhour);

%% (2) 初始化
a_pred = zeros(size(z));
a_pred(1) = z(1); % 第一个观测点作为起始

% 周期步长
if isfield(p, 'step_size')
    step = p.step_size;
else
    step = 1000;
end

%% (3) 初始化裂纹几何
% 使用第一次观测时权重最大粒子的几何信息作为起始点
if isfield(p, 'init_geometry') && ~isempty(p.init_geometry)
    % 使用保存的真实几何信息
    yRegSet = p.init_geometry.yRegSet;
    zRegSet = p.init_geometry.zRegSet;
    n_nodes = length(yRegSet);
    % fprintf('  [pred_a_N] 使用保存的初始几何信息 (n_nodes=%d, z_end=%.4f)\\n', ...
        % n_nodes, zRegSet(end));
else
    % 如果没有保存的几何信息，使用默认贯穿裂纹初始化
    n_nodes = 21;
    yRegSet = linspace(11, 9, n_nodes);
    zRegSet = repmat(z(1), 1, n_nodes);
    warning('pred_a_N: 未找到初始几何信息，使用默认贯穿裂纹初始化');
end

% 加载 POD 模型数据
load('pod_models.mat');

% 加载测试误差集
testErrSet = [0.019803420755871 0.058377623733943 0.013343080193008 ...
    0.009221345729625 0.037283991994545 0.011408674471790 ...
    0.082711915490345 0.032375406062069 0.041104885753232 ...
    0.057544490601307];

% 阶段映射
stage = [1, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 8];

% 分裂状态
SPLITTED = 0;

% 预加载分裂模型数据
Uinput_splitted_1 = Uinput_splitted{1};
averInput_splitted_1 = averInput_splitted{1};
Uinput_splitted_2 = Uinput_splitted{2};
averInput_splitted_2 = averInput_splitted{2};

%% (4) 迭代流程
N_now = N_check(1);
N_total_end = N_check(end);
next_check_idx = 2;

% 转换参数格式
% 注意: PARA0 = [k1, k2, log_theta1, theta3, theta2]
% 所以 f = [log_theta1, theta3, theta2]
log_theta1_ = f(1);  % log_theta1
theta3 = f(2);       % theta3 (不是 theta2!)
theta2 = f(3);       % theta2 (不是 theta3!)
k2 = k(2);

while N_now < N_total_end
    % 计算步长
    actual_step = min(step, N_check(next_check_idx) - N_now);
    
    %% 裂纹阶段判断
    a_up = zRegSet(end) - 30;
    a_down = zRegSet(1) - 30;
    
    m_index = 0;
    while m_index == 0
        m_index = getModelIndexFunc(a_up, a_down);
        if m_index == 0
            % 如果无法判断阶段，使用默认值
            m_index = 8; % 默认使用最后阶段
            break;
        end
    end
    
    %% 模型选择
    if SPLITTED && (m_index == 3 || m_index == 5)
        if m_index == 3
            curUinput = Uinput_splitted_1;
            curAverInput = averInput_splitted_1;
        elseif m_index == 5
            curUinput = Uinput_splitted_2;
            curAverInput = averInput_splitted_2;
        end
        m_name = sprintf('nn_stage%ds', stage(m_index));
    else
        curUinput = Uinput_integrated{stage(m_index)};
        curAverInput = averInput_integrated{stage(m_index)};
        m_name = sprintf('nn_stage%d', stage(m_index));
    end
    
    %% 提取载荷谱段
    idx_start = 2 * N_now + 1;
    idx_end = 2 * (N_now + actual_step);
    
    if idx_end > length(spectra)
        break;
    end
    
    loads_segment = spectra(idx_start : idx_end);
    
    %% 调用 a2aNew 更新裂纹状态
    try
        [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, ~] = ...
            a2aNew(yRegSet, zRegSet, loads_segment, m_name, ...
            curUinput, curAverInput, log_theta1_, theta2, theta3, k2, actual_step, testErrSet, 0);
    catch ME
        % 如果预测失败，返回无穷大的 MSE
        mse = inf;
        a_pred = nan(size(z));
        return;
    end
    
    %% 更新时钟
    N_now = N_now + actual_step;
    
    %% 记录检查点
    if N_now >= N_check(next_check_idx)
        a_pred(next_check_idx) = zRegSet(end); % 使用最后一个节点的 z 坐标
        next_check_idx = next_check_idx + 1;
        if next_check_idx > length(N_check)
            break;
        end
    end
end

%% (5) 计算 MSE
mse = mean((a_pred - z).^2) / 2;

end
