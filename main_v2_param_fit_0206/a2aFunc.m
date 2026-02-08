%% ===================================================================
%% 函数名称：a2aFunc
%% 功能描述：基于POD神经网络模型的裂纹扩展预测核心算法
%% ===================================================================
function [yRegSet, zRegSet, SPLITTED, log_theta1_, theta2, theta3, k2, deltaKSet] = ...
    a2aFunc(yRegSet, zRegSet, loads_segment, m_name, curUinput, curAverInput, log_theta1_, theta2, theta3, k2, step, testErrSet, particle_idx)

%A2AFUNC 基于当前时刻的裂纹状态和载荷谱段预测下一时刻的裂纹状态
%
% 输入参数：
%   loads_segment:    当前 step 内的载荷谱段 [Smin1, Smax1, Smin2, Smax2, ...]
% ... (其他参数同原函数)

%% ===================================================================
%% 参数初始化和预处理
%% ===================================================================
nRegPoint = length(yRegSet);   % 裂纹轮廓节点数量
yInputSet = yRegSet;
zInputSet = zRegSet;

%% POD降维投影
inputRegSet = [yRegSet, zRegSet]';
input = curUinput' * (inputRegSet - curAverInput);

%% ===================================================================
%% 循环级扩展计算 (与 pred_a_N 逻辑同步)
%% ===================================================================

% 1. 提取 Smin, Smax 序列
Smin = loads_segment(1:2:end);
Smax = loads_segment(2:2:end);
delta_sigma = Smax - Smin;

% 2. 计算基准 K 值 (使用神经网络，以 130MPa 为基准载荷)
% 注意：0.367 是原程序中的特定转换系数，此处保留以维持物理一致性
ref_stress = 130;
[K_ref] = sim_K_func(m_name, input, ref_stress * 0.367, testErrSet); 
K_ref_m = K_ref' / sqrt(1000); % 转换为 MPa*sqrt(m)

% 3. 逐循环计算增量并累加
da_total = zeros(1, nRegPoint);
bracket_all = zeros(length(Smin), nRegPoint);

theta1 = 10 ^ (log_theta1_);

% 计算每个循环的 K 序列
% deltaK = (delta_sigma / 130) * K_ref_m
% Kmax = (Smax / 130) * K_ref_m
for i = 1:length(Smin)
    dk_cycle = (delta_sigma(i) / ref_stress) .* K_ref_m;
    kmax_cycle = (Smax(i) / ref_stress) .* K_ref_m;
    
    R_ratio = 1 - dk_cycle ./ kmax_cycle;
    delta_K_th = k2 .* (1 - R_ratio);
    
    % 阈值检查
    idx_growth = (dk_cycle - delta_K_th) >= 0;
    
    % 计算该循环的增量 (m)
    bracket_term = max(kmax_cycle ./ k2 - 1, 0);
    da_cycle = theta1 .* (dk_cycle .^ theta2) .* (bracket_term .^ theta3);
    da_cycle(~idx_growth) = 0;
    
    da_total = da_total + da_cycle;
end

% 转换为 mm 并作为最终的 da
da = 1e3 * da_total;

% 为 POF 计算提供参考 ΔK
deltaKSet = mean(Smax / ref_stress) .* K_ref; 

%% ===================================================================
%% 特殊情况处理：小裂纹的各向同性扩展
%% ===================================================================
a_old = sqrt((yRegSet - 13).^2 + (zRegSet - 30).^2);
a_up = a_old(end); 

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

ksiRegSet = linspace(0, 1, nRegPoint);  % 节点参数化坐标 (0~1)
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

%% ===================================================================
%% 异常状态检查（鲁棒性增强）
%% ===================================================================
max_da = max(da);

%% 检查坐标点数量是否过少
if length(yNewSet) <= 1 || length(zNewSet) <= 1
    fprintf('  [坐标检查] 粒子%d裂纹坐标点数量过少（仅剩%d个点），max(da)=%.4f\n', ...
        particle_idx, min(length(yNewSet), length(zNewSet)), max_da);
    % 如果是极端参数导致的异常，标记为无效粒子供上层捕捉
    error('A2A:InvalidParticle', '粒子%d发现极端扩展(max_da=%.2f)导致点数不足，标记为无效粒子。', ...
        particle_idx, max_da);
    % if max_da > 3
    %     % 如果是极端参数导致的异常，标记为无效粒子供上层捕捉
    %     error('A2A:InvalidParticle', '粒子%d发现极端扩展(max_da=%.2f)导致点数不足，标记为无效粒子。', ...
    %         particle_idx, max_da);
    % else
    %     % 否则保留原有报错逻辑（通过 fprintf 显示详情并终止）
    %     fprintf('  剩余坐标点信息:\n');
    %     for k = 1:min(length(yNewSet), length(zNewSet))
    %         fprintf('    第%d个点: (z=%.6f, y=%.6f)\n', k, zNewSet(k), yNewSet(k));
    %     end
    %     % 绘制诊断图形
    %     plotCrackCoordinates(yInputSet, zInputSet, particle_idx, -1, -1);
    %     plotCrackCoordinates(yOrigSet, zOrigSet, particle_idx, -2, -2);
    %     plotCrackCoordinates(yNewSet, zNewSet, particle_idx, -3, -3);
    %     error('[坐标检查] 粒子%d裂纹坐标点数量不足（非极端扩展），无法继续计算！', particle_idx);
    % end
end

%% 检查边界违规
[hasViolation, violationInfo] = checkBoundaryViolation(yNewSet, zNewSet);
if hasViolation
    fprintf('  [边界检查] 粒子%d发现裂纹坐标超出边界，max(da)=%.4f\n', particle_idx, max_da);
    error('A2A:InvalidParticle', '粒子%d检测到极端扩展(max_da=%.2f)导致边界违规，标记为无效粒子。', ...
    particle_idx, max_da)
    % if max_da > 3
    %     % 如果是极端参数导致的异常，标记为无效粒子供上层捕捉
    %     error('A2A:InvalidParticle', '粒子%d检测到极端扩展(max_da=%.2f)导致边界违规，标记为无效粒子。', ...
    %         particle_idx, max_da);
    % else
    %     % 否则保留原有报错逻辑（通过 fprintf 显示详情并终止）
    %     fprintf('  违规坐标点详情:\n');
    %     for k = 1:length(violationInfo)
    %         fprintf('    第%d个点 (z=%.6f, y=%.6f) 在%s: %s\n', ...
    %             violationInfo(k).point, violationInfo(k).z, violationInfo(k).y, ...
    %             violationInfo(k).region, violationInfo(k).reason);
    %     end
    %     % 绘制诊断图形
    %     plotCrackCoordinates(yInputSet, zInputSet, particle_idx, -1, -1);
    %     plotCrackCoordinates(yOrigSet, zOrigSet, particle_idx, -2, -2);
    %     plotCrackCoordinates(yNewSet, zNewSet, particle_idx, -3, -3);
    %     error('[边界检查] 粒子%d检测到裂纹坐标超出物理边界（非极端扩展）！', particle_idx);
    % end
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

