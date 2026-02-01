function xparticles_new = resample_particles_independent(xparticles, weights_matrix)
% RESAMPLE_PARTICLES_INDEPENDENT 对每个状态维度独立进行重采样
% 
% 输入:
%   xparticles     - 所有的粒子状态矩阵 [N x D]
%   weights_matrix - 权重矩阵 [N x D]，每一列对应一个维度的似然权重
%
% 输出:
%   xparticles_new - 重采样后的粒子状态矩阵 [N x D]

    [N, D] = size(xparticles);
    xparticles_new = zeros(N, D);
    
    for d = 1:D
        % 归一化当前维度的权重
        w = weights_matrix(:, d);
        
        % 处理 NaN
        w(isnan(w)) = 0;
        
        sw = sum(w);
        if sw > 0 && ~isnan(sw) && ~isinf(sw)
            w = w / sw;
        else
            w = ones(N, 1) / N;
        end
        
        % 执行多项式重采样 (Multinomial Resampling)
        % 构造累积分布
        edges = [0; cumsum(w)];
        
        % 修复浮点数精度导致的非单调问题
        % 如果 edges(end) 略大于 1，手动强制设为 1 并确保前面不大于 1
        edges = min(edges, 1.0);
        % 确保单调不减
        for i = 2:length(edges)
            if edges(i) < edges(i-1)
                edges(i) = edges(i-1);
            end
        end
        
        % 均匀分布采样
        u = rand(N, 1);
        
        % 查找索引
        [~, idx] = histc(u, edges);
        
        % 安全处理索引：histc 可能返回 0 (未找到) 或 N+1 (正好等于最后一个边界)
        idx(idx == 0) = 1; 
        idx(idx > N) = N;
        
        % 组装当前维度
        xparticles_new(:, d) = xparticles(idx, d);
    end
end
