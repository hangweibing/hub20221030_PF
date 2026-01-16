function [outindex] = randomr(w)
% RANDOMR 粒子滤波重采样函数
% 输入: w - 粒子权重向量
% 输出: outindex - 重采样后的粒子索引

N = length(w);
outindex = zeros(1, N);

% 1. 计算累积权重（归一化确保总和为1）
cw = cumsum(w(:));
cw = cw / cw(end); 

% 2. 生成随机数（受全局 rng 种子控制）
% 注意：此处不使用 parfor 以保证结果的严格可重复性
% 且对于 N=1000 左右的规模，向量化或简单循环比并行更快
u = rand(1, N);

% 3. 寻找对应的索引
for i = 1:N
    idx = find(cw >= u(i), 1);
    if isempty(idx)
        outindex(i) = N;
    else
        outindex(i) = idx;
    end
end

end
