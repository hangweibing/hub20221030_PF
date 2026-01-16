%% ===================================================================
%% 函数名称：sim_K_func
%% 功能描述：基于神经网络模型预测应力强度因子范围
%% ===================================================================
function [deltaKSet] = sim_K_func(m_name, input, aver_delta_sigma, testErrSet)

% =========================================================================
% 输入参数：
%   m_name:           神经网络模型名称 (如 'nn_stage3', 'nn_stage5s')
%   input:            POD降维后的几何特征输入
%   aver_delta_sigma: 平均应力增量 (MPa)
%   testErrSet:       测试误差集 (10个误差值)
%
% 输出参数：
%   deltaKSet:        应力强度因子范围预测值 (考虑应力幅值的归一化结果)
%
% 功能说明：
%   通过动态调用预训练的神经网络模型，预测裂纹各节点的应力强度因子
%   并进行后处理以确保预测结果的合理性和稳定性
% =========================================================================

%% ===================================================================
%% 动态模型调用和预测
%% ===================================================================

% 构建动态函数调用命令
% 将模型名称转换为函数调用语句，如：KSet = nn_stage3(input)
commd = sprintf('KSet=%s(input);', m_name);

% 执行动态函数调用，获取神经网络预测结果
eval(commd);

%% ===================================================================
%% 误差参数选择
%% ===================================================================

% 从模型名称中提取阶段索引
% 例如：'nn_stage3' → '3', 'nn_stage5s' → '5s'
crackIndex = strrep(m_name, 'nn_stage', '');

% 根据阶段索引选择对应的测试误差
if length(crackIndex) == 1
    % 标准模型：直接使用数字索引
    testErr = testErrSet(str2num(crackIndex));
else
    % 分裂模型：使用特殊映射
    if crackIndex(1) == '3'
        testErr = testErrSet(9);  % 第3阶段分裂模型使用第9个误差
    elseif crackIndex(1) == '5'
        testErr = testErrSet(10); % 第5阶段分裂模型使用第10个误差
    end
end

%% ===================================================================
%% 预测结果后处理
%% ===================================================================

% 处理负值预测结果（物理上无意义）
% 将负的应力强度因子设置为20 MPa√m（保守估计）
KSet(KSet < 0) = 20;

% 使用线性插值方法移除异常值
% filloutliers函数自动检测并修复离群值
KSet = filloutliers(KSet, 'linear');

%% ===================================================================
%% 已注释的高级异常值处理方法
%% ===================================================================

% 详细的异常值检测和修复算法（已弃用）
% 该方法通过邻域插值修复异常预测值，但计算复杂度较高

% if ~isempty(outliers)
%     for k=1:length(outliers)
%         % 找到异常值在正常点群中的位置
%         if outliers(k)<min(normal_points)
%             neighbor_points=normal_points(find(normal_points>1,2));
%             KSet(outliers(k))=KSet(neighbor_points(1));
%         elseif outliers(k)>max(normal_points)
%             tmp_array=find(normal_points<num_K_points);
%             neighbor_point=normal_points(tmp_array(end));
%             KSet(outliers(k))=KSet(neighbor_point);
%         else
%             tmp_array=find(normal_points<outliers(k));
%             left_neighbor=normal_points(tmp_array(end));
%             right_neighbor=normal_points(find(normal_points>outliers(k),1));
%             neighbor_points=[left_neighbor,right_neighbor];
%             KSet(outliers(k))=interp1([KSet(neighbor_points(1)),KSet(neighbor_points(2))],[neighbor_points(1),neighbor_points(2)],outliers(k),'linear','extrap');
%         end
%     end
% end

%% ===================================================================
%% 简单阈值滤波方法（已弃用）
%% ===================================================================

% 基于平均值的简单异常值检测（已弃用）
% 该方法过于简单，可能误删有效预测结果

% for i=1:size(KSet,2)
%     if abs(KSet(i)-averK)>0.5*abs(averK)
%         if i==1
%             KSet(i)=KSet(2);
%         elseif i==size(KSet,2)
%             KSet(i)=KSet(size(KSet,2)-1);
%         else
%             KSet(i)=(KSet(i-1)+KSet(i+1))/2;
%         end
%     end
% end

%% ===================================================================
%% 调试信息输出（已注释）
%% ===================================================================

% 输出预测结果的统计信息（用于调试，已注释）
% fprintf('Kset mean: %f m_name: %s ',mean(KSet),m_name)

%% ===================================================================
%% 最终结果计算和归一化
%% ===================================================================

% 将预测的应力强度因子转换为实际的应力强度因子范围
% KSet是基于归一化应力(130MPa)训练的，需要根据实际应力幅值进行调整
deltaKSet = KSet .* aver_delta_sigma / 130;

end