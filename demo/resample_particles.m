function [xparticle_new, weight_new] = resample_particles(xparticle_curr, weights)
% RESAMPLE_PARTICLES Perform particle resampling
% xparticle_curr: Current particle population matrix
% weights:        Normalized particle weights

N = size(xparticle_curr, 1);

% 使用 randomr 函数进行系统重采样 (randomr.m 应该在路径中)
% 如果 randomr 不在 demo 文件夹，请确保 main 文件夹在 MATLAB 路径中
outindex = randomr(weights);

% 更新粒子群体
xparticle_new = xparticle_curr(outindex, :);

% 重置权重为均匀分布
weight_new = 1/N * ones(N, 1);
end
