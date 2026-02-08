function [outindex] = randomr(w, seed_offset)
% RANDOMR 按权重进行随机重采样，使用线程安全的随机数生成
% 输入参数：
%   w - 权重向量
%   seed_offset - 种子偏移量（可选），用于确保可重现性
%
% 如果提供了seed_offset，则使用确定性种子；否则使用全局随机数

N = length(w);

if nargin >= 2
    % 使用确定性种子确保可重现性
    stream = RandStream('mt19937ar', 'Seed', seed_offset);
    outindex = zeros(1, N);
    for i = 1:N
        outindex(i) = find(rand(stream) <= cumsum(w), 1);
    end
else
    % 使用全局随机数（保留向后兼容性）
    parfor i = 1:N
        outindex(i) = find(rand <= cumsum(w), 1);
    end
end
end

