%% ===================================================================
%% 函数名称：a2aNew
%% 功能描述：基于神经网络模型预测裂纹几何形状和材料参数的演化
%% ===================================================================
function [yRegSet, zRegSet, SPLITTED, logCstar, gamma] = ...
         a2aNew(yRegSet, zRegSet, aver_delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step, testErrSet)

%A2ANEW 基于当前时刻的裂纹状态和对应模型预测下一时刻的裂纹状态
%
% 输入参数：
%   yRegSet:         当前时刻y坐标集 (21个节点)
%   zRegSet:         当前时刻z坐标集 (21个节点)
%   aver_delta_sigma: 平均应力增量 (MPa)
%   m_name:          模型名称 (如 'nn_stage3', 'nn_stage5s' 等)
%   curUinput:       分裂模型的U输入数据
%   curAverInput:    分裂模型的平均输入数据
%   logCstar:        参数 logC* (当前值)
%   gamma:           参数 γ (当前值)
%   step:            时间步长
%   testErrSet:      测试误差集 (用于模型不确定性)
%
% 输出参数：
%   yRegSet:         预测的下一时刻y坐标集
%   zRegSet:         预测的下一时刻z坐标集
%   SPLITTED:        裂纹是否发生分裂 (0/1)
%   logCstar:        更新的Paris定律参数 logC*
%   gamma:           更新的Paris定律参数 γ
%
% 算法流程：
%   1. 基于当前裂纹几何形状和应力增量
%   2. 使用对应阶段的神经网络模型进行预测
%   3. 输出下一时刻的裂纹状态和材料参数
%   4. 同时判断是否发生裂纹分裂现象

%% 核心预测逻辑
% 调用a2aFunc函数执行实际的神经网络预测

%% 注释：早期版本的条件判断逻辑（已弃用）
% 原设计：当裂纹深度超过35mm时停止扩展
% if max(zRegSet) >= 35
%     % 保持当前状态不变，不再扩展
%     yRegSet = yRegSet;
%     zRegSet = zRegSet;
%     SPLITTED = 0;
%     logCstar = logCstar;
%     gamma = gamma;
% else
%     % 调用预测函数
%     [yRegSet, zRegSet, SPLITTED, logCstar, gamma] = ...
%         a2aFunc(yRegSet, zRegSet, delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step);
% end

%% 执行神经网络预测
[yRegSet, zRegSet, SPLITTED, logCstar, gamma] = ...
    a2aFunc(yRegSet, zRegSet, aver_delta_sigma, m_name, curUinput, curAverInput, logCstar, gamma, step, testErrSet);