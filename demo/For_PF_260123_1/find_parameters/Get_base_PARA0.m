function [selectedParams] = Get_base_PARA0( A, maxK)
    warning off;
    % 第 1 列是 loss，第 2–6 列是拟合出的 5 个参数
    loss   = A(:,1);          
    params = A(:,2:end);
    N      = size(params,1);  % 样本数量

    % 处理样本不足的情况
    if N <= 1
        selectedParams = params;
        return;
    end
    
    %% 1. 先把参数标准化（零均值、单位方差） 
    Z = (params - mean(params,1))./std(params,0,1);  
    Z(isnan(Z)) = 0;  % 处理标准差为0导致的NaN    

    %% 2. 动态调整聚类范围 —— 这里用 silhouette
    %maxK   = min(10, size(Z,1)-1);   % 不会超过 10，也不会超过 N-1
    kRange = 2:maxK;    % 至少分 2 类才有意义
 

    %% 3. 自动决定最佳聚类数（带错误处理）
    try
        eva = evalclusters(Z, 'kmeans', 'silhouette', 'KList', kRange);
        bestK = eva.OptimalK;
    catch
        % 聚类失败时使用距离法选点
        D = pdist2(Z, Z);                      % 距离矩阵
        [~, centerIdx] = min(sum(D, 2));       % 找中心点
        [~, farIdx] = max(D(centerIdx, :));    % 找离中心最远的点
        
        % 确保至少选2个点
        if N >= 2
            selectedIndices = [centerIdx; farIdx];
        else
            selectedIndices = centerIdx;
        end
        
        selectedParams = params(selectedIndices, :);
        return;
    end

    %% 4. 执行K-means聚类
    rng default  % 保证可复现
    [idx, ~] = kmeans(Z, bestK, 'Replicates', 20, 'Display', 'off');

    %% 5. 每类选代表点（优先loss最小，次选近中心点）
    % bestRows = zeros(bestK,1);
    % for k = 1:bestK
    %     clusterIndices = find(idx == k);
    % 
    %     % 优先选loss最小的点
    %     [minLoss, minIdx] = min(loss(clusterIndices));
    % 
    %     % 特殊情况处理：如果所有loss相同，选离中心最近的点
    %     if all(loss(clusterIndices) == minLoss)
    %         [~, minDistIdx] = min(pdist2(C(k,:), Z(clusterIndices,:)));
    %         bestRows(k) = clusterIndices(minDistIdx);
    %     else
    %         bestRows(k) = clusterIndices(minIdx);
    %     end
    % end

    bestRows = nan(bestK,1);
    for k = 1:bestK
        rowsInCluster   = find(idx == k);
        [~,localMinLoc] = min(loss(rowsInCluster));
        bestRows(k)     = rowsInCluster(localMinLoc);
    end


    %% 6. 输出代表行（含5参数）
    selectedParams = params(bestRows, :);

end

