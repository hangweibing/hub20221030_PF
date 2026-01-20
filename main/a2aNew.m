%% ===================================================================
%% 函数名称：a2aNew
%% 功能描述：a2aFunc的包装函数，增加deltaKSet输出用于POF计算
%% ===================================================================
function [yRegSet, zRegSet, SPLITTED, logD, A, delta_kthr, p, deltaKSet] = ...
         a2aNew(yRegSet, zRegSet, aver_delta_sigma, aver_R, m_name, curUinput, curAverInput, logD, A, delta_kthr, p, step, testErrSet, particle_idx)
% A2ANEW a2aFunc的包装函数，返回应力强度因子用于失效概率计算
%
% 该函数是 a2aFunc 的简单包装，确保返回 deltaKSet 用于后续的
% 失效概率(POF)计算。所有参数和功能与 a2aFunc 完全相同。
%
% 输入输出参数：与 a2aFunc 完全相同，参见 a2aFunc.m 文档
%
% 新增输出：
%   deltaKSet: 应力强度因子范围，用于POF计算

%% 直接调用 a2aFunc
[yRegSet, zRegSet, SPLITTED, logD, A, delta_kthr, p, deltaKSet] = ...
    a2aFunc(yRegSet, zRegSet, aver_delta_sigma, aver_R, m_name, ...
            curUinput, curAverInput, logD, A, delta_kthr, p, step, testErrSet, particle_idx);

end
