%% ===================================================================
%% 函数名称：a2aFunc
%% 功能描述：基于POD神经网络模型的裂纹扩展预测核心算法
%% ===================================================================
function [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, deltaKSet] = ...
         a2aFunc(yRegSet, zRegSet, aver_delta_sigma, aver_R, aver_Smax, m_name, curUinput, curAverInput, log_theta1_, theta2, theta3, k2, step, testErrSet, particle_idx)

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
%   aver_R:           平均应力比 R = Smin/Smax
%   aver_Smax:        平均最大应力 (MPa)
%   m_name:          模型名称 (如 'nn_stage3', 'nn_stage5s' 等)
%   curUinput:       分裂模型的U输入矩阵 (POD基向量)
%   curAverInput:    分裂模型的平均输入数据
%   log_theta1_:     新的模型参数 theta1 (对数形式)
%   theta2:          新的模型参数 theta2
%   theta3:          新的模型参数 theta3
%   k2:              新的模型参数 k2
%   step:            时间步长 (默认值为1)
%   testErrSet:      测试误差集 (用于模型不确定性)
%
% 输出参数：
%   yRegSet:         预测的下一时刻y坐标集
%   zRegSet:         预测的下一时刻z坐标集
%   SPLITTED:        裂纹是否发生分裂 (0/1)
%   log_theta1_:     更新的模型参数 theta1
%   theta2:          更新的模型参数 theta2
%   theta3:          更新的模型参数 theta3
%   k2:              更新的模型参数 k2
%   deltaKSet:       应力强度因子范围（用于POF计算）
%

%% ===================================================================
%% 参数初始化和预处理
%% ===================================================================

% 参数默认值处理
if nargin <= 11
    step = 1;  % 默认时间步长
end

%% 基本参数计算
nRegPoint = length(yRegSet);   % 裂纹轮廓节点数量

% 保存最原始的输入坐标，用于边界违规时的对比分析
yInputSet = yRegSet;
zInputSet = zRegSet;


%% POD降维投影
% 将裂纹几何形状投影到POD基空间进行降维
inputRegSet = [yRegSet, zRegSet]';  % 组合y和z坐标为输入矩阵
input = curUinput' * (inputRegSet - curAverInput);  % POD投影，每一列表示一个坐标

%% ===================================================================
%% 神经网络预测K（！！！！！注意单位！！！！！）
%% ===================================================================
% 使用神经网络计算 deltaK Kmax
[deltaKSet] = sim_K_func(m_name, input, aver_delta_sigma * 0.367, testErrSet);
deltaKSet_vec = deltaKSet';  
deltaKSet_m = deltaKSet_vec / sqrt(1000);  

[Kmax] = sim_K_func(m_name, input, aver_Smax * 0.367, testErrSet);
Kmax = Kmax / sqrt(1000);
Kmax = Kmax';  % 转换为行向量，与deltaKSet_m保持一致的维度

% deltaKSet_vec = deltaKSet';  % 转换为行向量
% deltaKSet_m = deltaKSet_vec / sqrt(1000);  % 从mm单位转换为m单位
% Kmax = deltaKSet_m ./ (1 - aver_R) ;  % Kmax = ΔK / (1 - R)
%% 裂纹几何参数计算
ksiRegSet = linspace(0, 1, nRegPoint);  % 节点参数化坐标 (0~1)
y_ini = 13;                             % 初始圆心y坐标
z_ini = 30;                             % 初始圆心z坐标
a_old = sqrt((yRegSet - y_ini).^2 + (zRegSet - z_ini).^2);  % 当前裂纹尺寸

%% 新的裂纹扩展模型计算
% da = step * theta1 * (delta)^theta2 * (Kmax/k2 - 1)^theta3
% 其中：da-裂纹扩展量，delta-应力增量，Kmax-最大应力强度因子
%      theta1, theta2, theta3, k2 - 模型参数

% 计算新的模型公式
% 注意：需要确保括号内的值不为负数
bracket_term = Kmax ./ k2 - 1;
bracket_term = max(bracket_term, 0);  % 确保不为负数

% 新的模型公式：将对数参数转换为实际参数值
theta1 = 10 ^ (log_theta1_);  % 从对数形式转换为实际参数值
da = step * theta1 .* (deltaKSet_m .^ theta2) .* (bracket_term .^ theta3);
da = 1e3 * da; % 单位转换

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
% 保存约束处理前的原始坐标，用于对比分析
yOrigSet = yNewSet;
zOrigSet = zNewSet;
% 根据边界条件调整裂纹形状，同时检测是否发生分裂
[yNewSet, zNewSet, SPLITTED] = addConstraintNewSatgeFunc(yNewSet, zNewSet);

%% 检查坐标点数量是否过少
if length(yNewSet) <= 1 || length(zNewSet) <= 1
    fprintf('  [坐标检查] 粒子%d裂纹坐标点数量过少（仅剩%d个点），正在绘制图像并终止程序...\n', ...
        particle_idx, min(length(yNewSet), length(zNewSet)));
    fprintf('  剩余坐标点信息:\n');
    for k = 1:min(length(yNewSet), length(zNewSet))
        fprintf('    第%d个点: (z=%.6f, y=%.6f)\n', k, zNewSet(k), yNewSet(k));
    end
    % 绘制最原始输入的图形
    plotCrackCoordinates(yInputSet, zInputSet, particle_idx, -1, -1);  % 包含粒子编号
    % 绘制处理前的原始图形
    plotCrackCoordinates(yOrigSet, zOrigSet, particle_idx, -2, -2);  % 包含粒子编号
    % 绘制约束处理后的图形
    plotCrackCoordinates(yNewSet, zNewSet, particle_idx, -3, -3);   % 包含粒子编号
    error('[坐标检查] 粒子%d裂纹坐标点数量不足，无法继续计算！', particle_idx);
end

%% 检查边界违规
[hasViolation, violationInfo] = checkBoundaryViolation(yNewSet, zNewSet);
if hasViolation
    fprintf('  [边界检查] 粒子%d发现裂纹坐标超出边界，正在绘制图像并终止程序...\n', particle_idx);
    fprintf('  违规坐标点详情:\n');
    for k = 1:length(violationInfo)
        fprintf('    第%d个点 (z=%.6f, y=%.6f) 在%s: %s\n', ...
            violationInfo(k).point, violationInfo(k).z, violationInfo(k).y, ...
            violationInfo(k).region, violationInfo(k).reason);
    end
    % 绘制最原始输入的图形
    plotCrackCoordinates(yInputSet, zInputSet, particle_idx, -1, -1);  % 包含粒子编号
    % 绘制处理前的原始图形
    plotCrackCoordinates(yOrigSet, zOrigSet, particle_idx, -2, -2);  % 包含粒子编号
    % 绘制约束处理后的图形
    plotCrackCoordinates(yNewSet, zNewSet, particle_idx, -3, -3);   % 包含粒子编号
    error('[边界检查] 粒子%d检测到裂纹坐标超出物理边界！程序终止以确保计算质量。', particle_idx);
end

%% 几何正则化
% t_regular = tic;  % 性能统计已注释
% 重新生成规则的裂纹轮廓，确保几何约束满足且没有自相交
[yRegSet, zRegSet, ~, SPLITTED] = crackRegular5Func(yNewSet, zNewSet, nRegPoint, 'false');



%% 输出参数（当前版本保持材料参数不变）
log_theta1_ = log_theta1_;  % 模型参数theta1 (保持不变)
theta2 = theta2;            % 模型参数theta2 (保持不变)
theta3 = theta3;            % 模型参数theta3 (保持不变)
k2 = k2;                    % 模型参数k2 (保持不变)
% deltaKSet 已在第62行计算，作为输出返回用于POF计算
end

