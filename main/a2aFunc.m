%% ===================================================================
%% 函数名称：a2aFunc
%% 功能描述：基于POD神经网络模型的裂纹扩展预测核心算法
%% ===================================================================
function [yRegSet, zRegSet, SPLITTED, logCstar, gamma] = ...
         a2aFunc(yRegSet, zRegSet, aver_delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step, testErrSet)

%A2AFUNC 基于当前时刻的裂纹状态和对应模型预测下一时刻的裂纹状态
%
% 算法原理：
%   采用POD（Proper Orthogonal Decomposition）降维和神经网络相结合的方法
%   预测裂纹在循环载荷作用下的几何形状演化
%
% 输入参数：
%   yRegSet:         当前时刻y坐标集 (21个节点)
%   zRegSet:         当前时刻z坐标集 (21个节点)
%   aver_delta_sigma: 平均应力增量 (MPa)
%   m_name:          模型名称 (如 'nn_stage3', 'nn_stage5s' 等)
%   curUinput:       分裂模型的U输入矩阵 (POD基向量)
%   curAverInput:    分裂模型的平均输入数据
%   logCstar:        Paris定律参数 logC* (当前值)
%   gamma:           Paris定律参数 γ (当前值)
%   step:            时间步长 (默认值为1)
%   testErrSet:      测试误差集 (用于模型不确定性)
%
% 输出参数：
%   yRegSet:         预测的下一时刻y坐标集
%   zRegSet:         预测的下一时刻z坐标集
%   SPLITTED:        裂纹是否发生分裂 (0/1)
%   logCstar:        更新的Paris定律参数 logC*
%   gamma:           更新的Paris定律参数 γ
%
% 核心技术：
%   1. POD (Proper Orthogonal Decomposition) 降维
%   2. 神经网络应力强度因子预测
%   3. Paris疲劳裂纹扩展定律
%   4. 几何约束和正则化处理
%% ===================================================================
%% 参数初始化和预处理
%% ===================================================================

% 参数默认值处理
if nargin <= 8
    step = 1;  % 默认时间步长
end

%% 基本参数计算
Cstar = 10^logCstar;           % Paris定律参数C* (从对数形式转换)
nRegPoint = length(yRegSet);   % 裂纹轮廓节点数量

%% ===================================================================
%% POD投影和神经网络预测
%% ===================================================================

%% POD降维投影
% t_pod = tic;  % 性能统计已注释
% 将裂纹几何形状投影到POD基空间进行降维
inputRegSet = [yRegSet, zRegSet]';  % 组合y和z坐标为输入矩阵
input = curUinput' * (inputRegSet - curAverInput);  % POD投影，每一列表示一个坐标
% t_pod_time = toc(t_pod);  % 性能统计已注释

%% 神经网络应力强度因子预测
% t_nn = tic;  % 性能统计已注释
% 使用训练好的神经网络模型预测各节点的应力强度因子范围
[deltaKSet] = sim_K_func(m_name, input, aver_delta_sigma, testErrSet);
% t_nn_time = toc(t_nn);  % 性能统计已注释
%% ===================================================================
%% Paris疲劳裂纹扩展定律计算
%% ===================================================================

% t_paris = tic;  % 性能统计已注释
%% 裂纹几何参数计算
ksiRegSet = linspace(0, 1, nRegPoint);  % 节点参数化坐标 (0~1)
y_ini = 13;                             % 初始圆心y坐标
z_ini = 30;                             % 初始圆心z坐标
a_old = sqrt((yRegSet - y_ini).^2 + (zRegSet - z_ini).^2);  % 当前裂纹尺寸

%% 随机扰动项（考虑材料和环境不确定性）
variance = 0.1;                         % 扰动方差
mu = -variance/2;                       % 均值调整（保持期望值为1）
omega = normrnd(mu, variance, 1, nRegPoint);  % 对数正态随机扰动

%% Paris定律裂纹扩展计算
% da/dN = C*(ΔK)^γ
% 其中：da-裂纹扩展量，dN-循环次数，ΔK-应力强度因子范围
da = step * Cstar .* a_old.^(1-gamma/2) .* deltaKSet'.^gamma .* exp(omega);
% t_paris_time = toc(t_paris);  % 性能统计已注释

%% ===================================================================
%% 特殊情况处理：小裂纹的各向同性扩展
%% ===================================================================

% t_small_crack = tic;  % 性能统计已注释
%% 判断是否为小裂纹阶段
a_up = a_old(end);  % 上表面裂纹尺寸

if(a_up <= 1.5)
    %% 小裂纹的圆形扩展模式
    % 对于小裂纹，假设近似圆形扩展，各方向扩展速率相近

    % 注释：早期版本使用均匀扩展
    % da1 = ones(1, nRegPoint) * da(1);   % 所有节点使用第一个节点的扩展量
    % da2 = ones(1, nRegPoint) * da(end); % 所有节点使用最后一个节点的扩展量

    % 当前版本：使用起点和终点的扩展量
    da1 = da(1);     % y方向基准扩展量
    da2 = da(end);   % z方向基准扩展量

    %% 计算圆形扩展的位移分量
    da_y_new = zeros(1, nRegPoint);  % y方向位移
    da_z_new = zeros(1, nRegPoint);  % z方向位移
    thetas = linspace(3/2*pi, 2*pi, nRegPoint);  % 角度分布（下半圆）

    for i = 1:nRegPoint
        % 将扩展量投影到圆形轨迹上
        da_y_new(i) = da1 * sin(thetas(i));  % y方向分量
        da_z_new(i) = da2 * cos(thetas(i));  % z方向分量
    end

    % 计算合位移（欧几里得距离）
    da_new = sqrt(da_y_new.^2 + da_z_new.^2);
    da = da_new;  % 更新扩展量
end
% t_small_crack_time = toc(t_small_crack);  % 性能统计已注释

%% ===================================================================
%% 几何形状更新
%% ===================================================================

%% 计算法向量（确定扩展方向）
% t_normal = tic;  % 性能统计已注释
normalNormalizeVectorSet = CalNormalVector(yRegSet, zRegSet, ksiRegSet);
% 计算每个节点的法向量，用于确定裂纹扩展的方向
% t_normal_time = toc(t_normal);  % 性能统计已注释

%% 计算坐标增量
% t_coord_update = tic;  % 性能统计已注释
yIncrSet = da .* normalNormalizeVectorSet(1, :);  % y方向增量
zIncrSet = da .* normalNormalizeVectorSet(2, :);  % z方向增量

%% 更新裂纹几何形状
yNewSet = yRegSet + yIncrSet;  % 新的y坐标集
zNewSet = zRegSet + zIncrSet;  % 新的z坐标集
% t_coord_update_time = toc(t_coord_update);  % 性能统计已注释

%% ===================================================================
%% 几何约束和正则化处理
%% ===================================================================

%% 边界约束处理
% t_constraint = tic;  % 性能统计已注释
% 根据边界条件调整裂纹形状，同时检测是否发生分裂
[yNewSet, zNewSet, SPLITTED] = addConstraintNewSatgeFunc(yNewSet, zNewSet);

% 检测复数并可视化（调试用）
if ~isreal(yNewSet) || ~isreal(zNewSet)
    warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
    plotCrackCoordinates(yNewSet, zNewSet, m_name);
end

% t_constraint_time = toc(t_constraint);  % 性能统计已注释

%% 几何正则化
% t_regular = tic;  % 性能统计已注释
% 重新生成规则的裂纹轮廓，确保几何约束满足且没有自相交
% 输入参数 'false' 表示不强制执行某些特殊约束
[yRegSet, zRegSet, ~, SPLITTED] = crackRegular5Func(yNewSet, zNewSet, nRegPoint, 'false');
% t_regular_time = toc(t_regular);  % 性能统计已注释


%% 输出参数（当前版本保持材料参数不变）
logCstar = logCstar;  % Paris定律参数logC* (保持不变)
gamma = gamma;        % Paris定律参数γ (保持不变)
end

