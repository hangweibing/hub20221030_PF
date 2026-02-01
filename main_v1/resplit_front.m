%% ===================================================================
%% 函数名称：resplit_front
%% 功能描述：处理裂纹前缘与下圆角的相交情况
%% 输入参数：
%%   yNewSet    - 裂纹的y坐标数组
%%   zNewSet    - 裂纹的z坐标数组
%%   downCenter - 下圆角圆心坐标 [x, z]
%%   r          - 下圆角半径
%%
%% 输出参数：
%%   yNewSet    - 处理后的y坐标数组（可能包含圆弧交点）
%%   zNewSet    - 处理后的z坐标数组（可能包含圆弧交点）
%%   splitted   - 是否发生裂纹分裂的标志（true/false）
%% ===================================================================

function [yNewSet, zNewSet, splitted] = resplit_front(yNewSet, zNewSet, downCenter, r)

%% 初始化参数
splitted = false;                    % 初始化分裂标志
n_reg_point = length(yNewSet);       % 获取节点总数
judge = zeros(1, n_reg_point);       % 初始化判断数组

%% 检查节点与圆角的相交情况
for i = 1:length(yNewSet)
    % 判断第i个节点是否在圆角区域内
    judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter, r);
end

%% 查找与圆角相交的节点索引
intersected_coords = find(judge == true);

%% 处理相交情况
if (~isempty(intersected_coords))
    % 存在与圆角相交的节点
    splitted = true;

    % 获取相交区域的边界索引
    left_coord = intersected_coords(1);     % 最左侧相交点索引
    right_coord = intersected_coords(end);  % 最右侧相交点索引

    %% 根据左侧边界位置进行不同处理
    if left_coord ~= 1
        %% 情况1：左侧边界不在起点，进行双侧边界点计算
        % 采用 projection 方法获取边界点，比求解直线圆弧交点更稳定
        [leftInsecY, leftInsecZ] = getEdgeYbyZFunc(zNewSet(left_coord), 'down', yNewSet(left_coord));

        % 采用 projection 方法获取边界点
        [rightInsecY, rightInsecZ] = getEdgeYbyZFunc(zNewSet(right_coord), 'down', yNewSet(right_coord));

        % 在左右交点之间生成10个插值点
        y_add = linspace(leftInsecY, rightInsecY, 10);
        z_add = linspace(leftInsecZ, rightInsecZ, 10);

        % 将这些插值点投影到下圆弧上
        for i = 1:length(z_add)
            [y_add(i), z_add(i)] = getEdgeYbyZFunc(z_add(i), 'down', y_add(i));
        end

        % 重新组合裂纹坐标：左侧段 + 圆弧段 + 右侧段
        yNewSet = [yNewSet(1:left_coord-1), y_add, yNewSet(right_coord+1:end)];
        zNewSet = [zNewSet(1:left_coord-1), z_add, zNewSet(right_coord+1:end)];

    else
        %% 情况2：左侧边界在起点，只处理右侧交点
        % 当下断点已经不在下边界时，只取上半段，需要对节点进行修正
        % 但下端点不一定在下圆弧上，有可能已经到达右下边界

        if zNewSet(right_coord) <= 37.82842712
            % 右侧点仍在圆弧范围内，采用 projection 方法获取边界点
            [rightInsecY, rightInsecZ] = getEdgeYbyZFunc(zNewSet(right_coord), 'down', yNewSet(right_coord));
        else
            % 右侧点已超出圆弧范围，进入直线边界区域 (y=9)
            % 同样采用 projection 方法获取边界点，确保与圆弧段处理逻辑一致并保持稳定性
            [rightInsecY, rightInsecZ] = getEdgeYbyZFunc(zNewSet(right_coord), 'down');
        end

        % 重置分裂标志（只处理一侧边界）
        splitted = false;

        % 重新组合裂纹坐标：交点 + 右侧段
        yNewSet = [rightInsecY, yNewSet(right_coord+1:end)];
        zNewSet = [rightInsecZ, zNewSet(right_coord+1:end)];
    end
end

%% 函数结束
end
