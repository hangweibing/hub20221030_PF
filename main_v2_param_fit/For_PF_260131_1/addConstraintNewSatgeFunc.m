%% ===================================================================
%% 函数名称：addConstraintNewSatgeFunc
%% 功能描述：根据工件几何约束调整裂纹轮廓形状
%% ===================================================================
function [yNewSet, zNewSet, splitted] = addConstraintNewSatgeFunc(yNewSet, zNewSet)

% 保存原始数据用于调试
origy = yNewSet;
origz = zNewSet;
if ~isreal(origy) || ~isreal(origz)
    warning('检测到复数坐标：origy 或 origz 包含复数');
    origy
    origz
end

% =========================================================================
% 输入参数：
%   yNewSet:         预测的裂纹y坐标集 (21个节点)
%   zNewSet:         预测的裂纹z坐标集 (21个节点)
%
% 输出参数：
%   yNewSet:         调整后的裂纹y坐标集
%   zNewSet:         调整后的裂纹z坐标集
%   splitted:        裂纹是否发生分裂 (0/1)
%
% 功能说明：
%   根据工件几何边界约束调整裂纹轮廓，确保裂纹形状符合物理约束
%   处理裂纹与圆角边界相交的情况，通过几何重构保证连续性
% =========================================================================

%% 初始化参数
splitted = false;                    % 分裂标志，初始为false
n_reg_point = length(yNewSet);       % 节点总数
minYloca = find(yNewSet < 7);        % 查找y坐标小于7的位置（接近左边界）

%% 几何约束参数定义
% 下圆角几何参数 (靠近孔的下方)
downCenter1 = [6, 37.82842712];      % 下圆角圆心坐标
upCenter1 = [14, 37.82842712];       % 上圆角圆心坐标

% 其他几何约束点 (备用)
downCenter2 = [1, 58];               % 备用下圆角
upCenter2 = [19, 58];                % 备用上圆角

r = 3;                               % 圆角半径 (mm)

%% 已注释的调试代码
% [~,index]=sort(yNewSet);
% if(index~=[1:11])
%     index=[1:11];
% end
% [minY,minYloca]=min(yNewSet);
% figure
% plot(zNewSet,yNewSet,'r'),axis equal
%% ===================================================================
%% 几何约束处理分支 - 根据裂纹位置和形状应用不同的约束规则
%% ===================================================================

%% 分支1：
if(isempty(minYloca) && zNewSet(end) <= 35)
    % (1) 处理左边界 (z=30)
    min_z_loca = find(zNewSet < 30);
    if ~isempty(min_z_loca)
        % 过滤掉所有违规点
        last_violation = min_z_loca(end);

        % 直接投影 (增加一个点而不是替换)
        yEdgeStart = yNewSet(last_violation+1);
        yNewSet = [yEdgeStart, yNewSet(last_violation+1:end)];
        zNewSet = [30, zNewSet(last_violation+1:end)];

    elseif ~isempty(zNewSet) && zNewSet(1) > 30.0
        % 如果起始点不在左边界上，增加交点 (直接投影)
        yNewSet = [yNewSet(1), yNewSet];
        zNewSet = [30, zNewSet];
    else
        zNewSet(1) = 30;
    end

    % (2) 处理上边界 (y=13)
    max_y_loca = find(yNewSet > 13);
    if ~isempty(max_y_loca)
        % 过滤掉违规点
        first_violation = max_y_loca(1);
        if first_violation > 1
            % 直接投影 (增加一个点而不是替换)
            zEdgeEnd = zNewSet(first_violation-1);
            zNewSet = [zNewSet(1:first_violation-1), zEdgeEnd];
            yNewSet = [yNewSet(1:first_violation-1), 13];
        end
    elseif ~isempty(yNewSet) && yNewSet(end) < 12.999
        % 如果末端点不在上边界上，增加交点 (直接投影)
        zNewSet = [zNewSet, zNewSet(end)];
        yNewSet = [yNewSet, 13];
    else
        yNewSet(end) = 13;
    end


    % 新增处理逻辑：判断是否有点的z坐标超过35mm（由于预测步长较大，部分节点可能已经进入圆弧区域）
    over35_indices = find(zNewSet > 35);
    if ~isempty(over35_indices)
        % 选取最后一个超过35的点
        lastZidx = over35_indices(end);
        if lastZidx > 1
            zEnd = zNewSet(lastZidx);
            % 画一条垂直线，计算该z坐标对应的上边界y值 (垂直线与圆弧的交点)
            yEnd = getEdgeYbyZFunc(zEnd, 'up');

            % 将此交点作为新的末端，填充到数组中并截断
            zNewSet = [zNewSet(1:lastZidx-1), zEnd];
            yNewSet = [yNewSet(1:lastZidx-1), yEnd];
        end
    end


    %% 分支2：
elseif(~isempty(minYloca) && zNewSet(end) <= 35)

    % (1) 处理下边界 (y=7)
    min_y_loca = find(yNewSet < 7);
    if ~isempty(min_y_loca)
        % 过滤掉违规点
        last_violation = min_y_loca(end);
        if last_violation < length(yNewSet)
            % 计算交点 (增加一个基于插值的起始点)
            % zEdgeStart = zNewSet(last_violation+1); % 原来的直接投影
            % yNewSet = [7, yNewSet(last_violation+1:end)];
            % zNewSet = [zEdgeStart, zNewSet(last_violation+1:end)];
            zEdgeStart = interp1(yNewSet(last_violation:last_violation+1), zNewSet(last_violation:last_violation+1), 7, 'linear', 'extrap');
            yNewSet = [7, yNewSet(last_violation+1:end)];
            zNewSet = [zEdgeStart, zNewSet(last_violation+1:end)];
        end
    elseif ~isempty(yNewSet) && yNewSet(1) > 7.001
        % 如果起始点不在下边界上，增加交点 (直接投影)
        zNewSet = [zNewSet(1), zNewSet];
        yNewSet = [7, yNewSet];
    else
        yNewSet(1) = 7;
    end

    % (2) 处理上边界 (y=13)
    max_y_loca = find(yNewSet > 13);
    if ~isempty(max_y_loca)
        % 过滤掉违规点
        first_violation = max_y_loca(1);
        if first_violation > 1
            % 直接投影 (增加一个点而不是替换)
            zEdgeEnd = zNewSet(first_violation-1);
            zNewSet = [zNewSet(1:first_violation-1), zEdgeEnd];
            yNewSet = [yNewSet(1:first_violation-1), 13];
        end
    elseif ~isempty(yNewSet) && yNewSet(end) < 12.999
        % 如果末端点不在上边界上，增加交点 (直接投影)
        zNewSet = [zNewSet, zNewSet(end)];
        yNewSet = [yNewSet, 13];
    else
        yNewSet(end) = 13;
    end

    % 同理，处理分支2中可能过线的点
    over35_indices = find(zNewSet > 35);
    if ~isempty(over35_indices)
        lastZidx = over35_indices(end);
        if lastZidx > 1
            zEnd = zNewSet(lastZidx);
            % 垂直投影求交点
            yEnd = getEdgeYbyZFunc(zEnd, 'up');

            zNewSet = [zNewSet(1:lastZidx-1), zEnd];
            yNewSet = [yNewSet(1:lastZidx-1), yEnd];
        end
    end

    %% 分支3：
elseif(isempty(minYloca) && (35 < zNewSet(end) && zNewSet(end) < 37.82842712) && zNewSet(1) < 35)

    % 新方法：参考前面分支的过滤和边界计算方法
    % 对超出左边界的节点进行过滤
    min_z_loca = find(zNewSet < 30);  % 查找超过左侧边界(z<30)的点
    if ~isempty(min_z_loca)
        % 保留从最后一个小z值位置到末尾的部分
        zNewSet = zNewSet(min_z_loca(end)+1:end);
        yNewSet = yNewSet(min_z_loca(end)+1:end);
    end

    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 检查每个点是否超出圆弧边界
    violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35  % 注意这里的 z 范围

            % 计算该z坐标对应的上边界和下边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);
            % 检查点是否超出边界
            if y > y_upper_boundary || y < y_lower_boundary
                violation_indices = [violation_indices, i];
            end

        end
    end

    % 如果有违规点，保留到第一个违规点之前的数据
    if ~isempty(violation_indices)
        first_violation = violation_indices(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % 增加与上圆弧的交点 (直接使用圆弧投影逻辑，避免因z(end)<35导致被判定为直线段)
    dz_up = zNewSet(end) - upCenter(1);
    dy_up = yNewSet(end) - upCenter(2);
    dist_up = sqrt(dz_up^2 + dy_up^2);
    yEdgeEnd = upCenter(2) + radius * dy_up / (dist_up + eps);
    zEdgeEnd = upCenter(1) + radius * dz_up / (dist_up + eps);
    yNewSet = [yNewSet, yEdgeEnd];
    zNewSet = [zNewSet, zEdgeEnd];

    % 确保起始点在左边界 (增加一个交点而不是替换)
    if ~isempty(zNewSet) && zNewSet(1) > 30.0
        yNewSet = [yNewSet(1), yNewSet];
        zNewSet = [30, zNewSet];
    else
        zNewSet(1) = 30;
    end





    %% 分支4：
elseif (~isempty(minYloca) || ((zNewSet(1)>30) && (zNewSet(1)<35))) && (zNewSet(end)>35) && (zNewSet(end)<=37.82842712)

    % 新方法：按步骤进行边界处理
    % （1）当minYloca非空，需要过滤超出下边界的节点
    if ~isempty(minYloca)
        min_y_loca = find(yNewSet < 7);  % 查找超过下边界(y<7)的点
        if ~isempty(min_y_loca)
            % 保留从最后一个违规点之后的部分，跳过所有违规点
            last_violation = min_y_loca(end);
            if last_violation < length(yNewSet)
                % 计算交点 (在last_violation+1之前增加一个和下边界的交点)
                % zNewSet = zNewSet(last_violation+1:end); % 原来的直接截断方式
                % yNewSet = yNewSet(last_violation+1:end);
                % yNewSet(1) = 7; % 原来的直接设置方式
                zEdgeStart = interp1(yNewSet(last_violation:last_violation+1), zNewSet(last_violation:last_violation+1), 7, 'linear', 'extrap');
                zNewSet = [zEdgeStart, zNewSet(last_violation+1:end)];
                yNewSet = [7, yNewSet(last_violation+1:end)];
            else
                yNewSet(1) = 7;
            end
        end
    end

    if isempty(minYloca)
        yNewSet(1) = 7;  % 直接设置第一个点的y坐标为7
    end

    max_y_loca = find(zNewSet > 37.82842712 & yNewSet > 11);  % 查找z>37.82842712且y>11的点
    if ~isempty(max_y_loca)
        % 保留从开始到第一个违规点之前的数据
        first_violation = max_y_loca(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % 过滤超出圆弧边界的点（参考分支三）
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 检查每个点是否超出圆弧边界
    violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35 && z <= 37.82842712
            % 计算该z坐标对应的上边界和下边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);

            % 检查点是否超出边界
            if y > y_upper_boundary || y < y_lower_boundary
                violation_indices = [violation_indices, i];
            end
        end
    end

    % 如果有违规点，保留到第一个违规点之前的数据
    if ~isempty(violation_indices)
        first_violation = violation_indices(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % 设置末端点到上圆弧 (增加一个交点而不是替换)
    if ~isempty(zNewSet)
        [yEdgeEnd, zEdgeEnd] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));
        yNewSet = [yNewSet, yEdgeEnd];
        zNewSet = [zNewSet, zEdgeEnd];
    end

    % 检查是否与下圆角相交 - 简化处理方式（参考分支7）
    % 对每个点判断是否在下圆弧内
    n_points = length(yNewSet);
    in_circle = false(n_points, 1);
    for i = 1:n_points
        in_circle(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
    end

    % 检测是否出现圆内点，只要有任何一个点在圆弧内就触发split处理
    has_split_pattern = any(in_circle);

    if has_split_pattern
        % 出现了split，过滤掉最后一个在圆内的点及其之前的点
        last_true_idx = find(in_circle, 1, 'last');
        yNewSet = yNewSet(last_true_idx+1:end);
        zNewSet = zNewSet(last_true_idx+1:end);
        n_points = length(yNewSet);

        % 根据第一个剩余点的z坐标找到边界交点（由于分支4的限制，一定在圆弧内）
        if n_points > 0
            % 使用getEdgeYbyZFunc函数计算下边界交点 (增加一个交点而不是替换)
            [yEdgeStart, zEdgeStart] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));
            yNewSet = [yEdgeStart, yNewSet];
            zNewSet = [zEdgeStart, zNewSet];
        end
    end


    %% 分支5：
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712 && zNewSet(1) > 35 && zNewSet(1) <= 37.82842712)

    % 新方法：分别过滤下圆弧和上圆弧边界外的点
    % 圆心坐标和半径（与checkBoundaryViolation.m保持一致）
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % （1）过滤超出下圆弧边界的点
    lower_violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35 && z < 37.82842712
            % 计算该z坐标对应的下边界y值
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);

            % 检查点是否超出下边界
            if y < y_lower_boundary
                lower_violation_indices = [lower_violation_indices, i];
            end
        end
    end

    % 如果有下边界违规点，保留最后一个违规点之后的数据
    if ~isempty(lower_violation_indices)
        last_violation = lower_violation_indices(end);
        if last_violation < length(zNewSet)
            zNewSet = zNewSet(last_violation+1:end);
            yNewSet = yNewSet(last_violation+1:end);
        end
    end

    % （2）过滤超出上圆弧边界的点
    upper_violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35 && z < 37.82842712
            % 计算该z坐标对应的上边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);

            % 检查点是否超出上边界
            if y > y_upper_boundary
                upper_violation_indices = [upper_violation_indices, i];
            end
        end
    end

    % 如果有上边界违规点，保留到第一个违规点之前的数据
    if ~isempty(upper_violation_indices)
        first_violation = upper_violation_indices(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % （3）计算和上下圆弧的交点：起始点与下圆弧相交，末端点与上圆弧相交
    if ~isempty(zNewSet)
        % 起始点与下圆弧相交 (增加一个交点而不是替换)
        if length(zNewSet) > 1
            [yEdgeStart, zEdgeStart] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));
            yNewSet = [yEdgeStart, yNewSet];
            zNewSet = [zEdgeStart, zNewSet];
        end

        % 末端点与上圆弧相交 (增加一个交点而不是替换)
        [yEdgeEnd, zEdgeEnd] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));
        yNewSet = [yNewSet, yEdgeEnd];
        zNewSet = [zNewSet, zEdgeEnd];
    end

    %% 分支6：
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712) && (zNewSet(1) >= 37.82842712)

    min_y_loca = find(yNewSet < 9);
    if ~isempty(min_y_loca)
        % 保留从最后一个越过下边界点之后的部分
        last_boundary_cross = min_y_loca(end);
        if last_boundary_cross < length(yNewSet)
            zNewSet = zNewSet(last_boundary_cross+1:end);
            yNewSet = yNewSet(last_boundary_cross+1:end);
        end
    end

    % （2）过滤超出圆弧范围的点，保留第一个超出范围之前的点
    % 圆心坐标和半径（与checkBoundaryViolation.m保持一致）
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 检查圆弧区域内的点是否超出边界
    violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35 && z < 37.82842712
            % 计算该z坐标对应的上边界和下边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);

            % 检查点是否超出边界
            if y > y_upper_boundary || y < y_lower_boundary
                violation_indices = [violation_indices, i];
            end
        end
    end

    % 如果有违规点，保留到第一个违规点之前的数据
    if ~isempty(violation_indices)
        first_violation = violation_indices(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % （3）定位和下边界以及上圆弧的交点
    % 设置起始点到下边界（y=9） (增加一个交点而不是替换)
    if ~isempty(yNewSet)
        if yNewSet(1) > 9.0
            yNewSet = [9, yNewSet];
            zNewSet = [zNewSet(1), zNewSet];
        else
            yNewSet(1) = 9;
        end
    end

    % 设置末端点到上圆弧 (增加一个交点而不是替换)
    if ~isempty(zNewSet) && zNewSet(end) < 37.82842712
        [yEdgeEnd, zEdgeEnd] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));
        yNewSet = [yNewSet, yEdgeEnd];
        zNewSet = [zNewSet, zEdgeEnd];
    end

    %% 分支7：
elseif(zNewSet(end) > 37.82842712) && (zNewSet(1) <= 35)
    % 处理起始点调整
    if ~isempty(minYloca)
        % 过滤超过下边界的点，保留最后一个违规点之后的数据
        min_y_loca = find(yNewSet < 7);  % 查找超过下边界(y<7)的点
        if ~isempty(min_y_loca)
            % 保留从最后一个违规点之后的部分
            last_violation = min_y_loca(end);
            if last_violation < length(yNewSet)
                zNewSet = zNewSet(last_violation+1:end);
                yNewSet = yNewSet(last_violation+1:end);
            end
        end
        % 设置起始点到下边界
        if ~isempty(yNewSet)
            yNewSet = [7, yNewSet];
            zNewSet = [zNewSet(1), zNewSet];
        end
    else
        % 过滤超过左边界的点，保留最后一个违规点之后的数据
        min_z_loca = find(zNewSet < 30);  % 查找超过左侧边界(z<30)的点
        if ~isempty(min_z_loca)
            % 保留从最后一个违规点之后的部分
            last_violation = min_z_loca(end);
            if last_violation < length(zNewSet)
                zNewSet = zNewSet(last_violation+1:end);
                yNewSet = yNewSet(last_violation+1:end);
            end
        end
        % 设置起始点到左边界 (增加一个交点而不是替换)
        if ~isempty(zNewSet) && zNewSet(1) > 30.0
            yNewSet = [yNewSet(1), yNewSet];
            zNewSet = [30, zNewSet];
        else
            zNewSet(1) = 30;
        end
    end

    % 检查是否有点在上圆弧区域内
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标 [z, y]
    radius = 3;                         % 圆弧半径

    % 检查是否有点在上圆弧区域内（z在圆心z坐标附近，考虑半径范围）
    points_in_upper_arc = false;
    for i = 1:length(yNewSet)
        point = [zNewSet(i), yNewSet(i)];  % [z, y] 坐标
        if isPointInCircleFunc(point, upCenter, radius)
            points_in_upper_arc = true;
            break;
        end
    end

    if points_in_upper_arc

        % 查找超出上圆弧的点
        violation_indices = [];
        for i = 1:length(yNewSet)
            z = zNewSet(i);
            y = yNewSet(i);
            point = [z, y];

            % 如果点在圆弧区域但超出上边界
            if isPointInCircleFunc(point, upCenter, radius)
                % 计算该z坐标对应的上边界y值
                if abs(z - upCenter(1)) <= radius
                    y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
                    if y > y_upper_boundary
                        violation_indices = [violation_indices, i];
                    end
                end
            end
        end

        % 过滤违规点
        if ~isempty(violation_indices)
            first_violation = violation_indices(1);
            if first_violation > 1
                zNewSet = zNewSet(1:first_violation-1);
                yNewSet = yNewSet(1:first_violation-1);
            end
        end
    else

        max_y_loca = find(zNewSet >= 37.82842712 & yNewSet > 11);  % 查找z>=37.82842712且y>11的点
        if ~isempty(max_y_loca)
            % 保留从开始到第一个违规点之前的数据
            first_violation = max_y_loca(1);
            if first_violation > 1
                zNewSet = zNewSet(1:first_violation-1);
                yNewSet = yNewSet(1:first_violation-1);
            end
        end
    end

    % 直接设置到上边界（y=11） (增加一个交点而不是替换)
    if ~isempty(yNewSet)
        if yNewSet(end) < 11.0
            yNewSet = [yNewSet, 11];
            zNewSet = [zNewSet, zNewSet(end)];
        else
            yNewSet(end) = 11;
        end
    end


    % 对每个点判断是否在该下边界违规区域 (圆弧内或超出直线边界)
    n_points = length(yNewSet);
    in_circle = false(n_points, 1);
    for i = 1:n_points
        % 检查下圆弧内点
        in_circle(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        % 检查 z > 37.82842712 时超出 y=9 边界的点
        if zNewSet(i) > 37.82842712 && yNewSet(i) < 9
            in_circle(i) = true;
        end
    end

    % 检测是否出现违规点，触发split处理
    has_split_pattern = any(in_circle);

    if has_split_pattern
        % 出现了split，过滤掉最后一个违规点及其之前的点
        last_true_idx = find(in_circle, 1, 'last');
        yNewSet = yNewSet(last_true_idx+1:end);
        zNewSet = zNewSet(last_true_idx+1:end);
        % 根据第一个剩余点的z坐标找到边界交点 (增加一个交点而不是替换)
        if ~isempty(yNewSet)
            [yEdgeStart, zEdgeStart] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));
            yNewSet = [yEdgeStart, yNewSet];
            zNewSet = [zEdgeStart, zNewSet];
        end
    end



    % 已注释的备用处理方法
    % [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)]);



    %% 分支8：
elseif(zNewSet(end) > 37.82842712 && zNewSet(1) > 35 && zNewSet(1) <= 37.82842712)
    % 新方法：系统化的边界过滤和交点计算
    % （1）过滤在圆弧外的点
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 检查圆弧区域内的点是否超出边界
    violation_indices = [];
    for i = 1:length(yNewSet)
        z = zNewSet(i);
        y = yNewSet(i);
        if z > 35 && z <= 37.82842712
            % 计算该z坐标对应的上边界和下边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);

            % 检查点是否超出边界
            if y > y_upper_boundary || y < y_lower_boundary
                violation_indices = [violation_indices, i];
            end
        end
    end

    % 如果有违规点，保留最后一个违规点之后的数据
    if ~isempty(violation_indices)
        last_violation = violation_indices(end);
        if last_violation < length(zNewSet)
            zNewSet = zNewSet(last_violation+1:end);
            yNewSet = yNewSet(last_violation+1:end);
        end
    end

    % 判断在z >= 37.82842712的区域，是否存在点在下边界（y=9）之外
    min_y_loca_region3 = find(zNewSet >= 37.82842712 & yNewSet < 9);  % 查找z>=37.82842712且y<9的点
    if ~isempty(min_y_loca_region3)
        % 保留最后一个违规点之后的数据
        last_violation_region3 = min_y_loca_region3(end);
        if last_violation_region3 < length(zNewSet)
            zNewSet = zNewSet(last_violation_region3+1:end);
            yNewSet = yNewSet(last_violation_region3+1:end);
        end
    end

    % 过滤上边界之外的点
    % 对于z >= 37.82842712的区域，上边界是y=11
    max_y_loca = find(zNewSet >= 37.82842712 & yNewSet > 11);  % 查找z>=37.82842712且y>11的点
    if ~isempty(max_y_loca)
        % 保留从开始到第一个违规点之前的数据
        first_violation = max_y_loca(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % 计算和上圆弧或上边界的交点
    if ~isempty(zNewSet)
        % 起始点：根据z坐标位置决定与哪个边界求交点
        [yEdgeStart, zEdgeStart] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));
        yNewSet = [yEdgeStart, yNewSet];
        zNewSet = [zEdgeStart, zNewSet];

        % 末端点：根据z坐标位置决定与哪个边界求交点
        if length(zNewSet) > 1
            [yEdgeEnd, zEdgeEnd] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));
            yNewSet = [yNewSet, yEdgeEnd];
            zNewSet = [zNewSet, zEdgeEnd];
        end
    end

    %% 分支9：
elseif zNewSet(end) > 37.82842712 && zNewSet(1) > 37.82842712  % 20220915 修改：上边界区域的处理
    % 原方法：清理边界并调整端点
    % min_y_loca = find(yNewSet < 9);  % 查找接近左边界(y<9)的点
    % if ~isempty(min_y_loca)
    %     zNewSet = zNewSet(min_y_loca(end):end);
    %     yNewSet = yNewSet(min_y_loca(end):end);
    % end
    %
    % max_y_loca = find(yNewSet > 11);  % 查找超过右边界(y>11)的点
    % if ~isempty(max_y_loca)
    %     zNewSet = zNewSet(1:max_y_loca(1));
    %     yNewSet = yNewSet(1:max_y_loca(1));
    % end
    %
    % [yNewSet(1), zNewSet(1)] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));
    % [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    % 新方法：参考前面的系统化边界过滤
    % 过滤超过下边界的点，保留最后一个超过下边界之后的点
    min_y_loca = find(yNewSet < 9);  % 查找超过下边界(y<9)的点
    if ~isempty(min_y_loca)
        % 保留从最后一个违规点之后的部分
        last_violation = min_y_loca(end);
        if last_violation < length(yNewSet)
            zNewSet = zNewSet(last_violation+1:end);
            yNewSet = yNewSet(last_violation+1:end);
        end
    end

    % 过滤超过上边界的点，保留第一个超过上边界之前的点
    max_y_loca = find(yNewSet > 11);  % 查找超过上边界(y>11)的点
    if ~isempty(max_y_loca)
        % 保留从开始到第一个违规点之前的数据
        first_violation = max_y_loca(1);
        if first_violation > 1
            zNewSet = zNewSet(1:first_violation-1);
            yNewSet = yNewSet(1:first_violation-1);
        end
    end

    % 计算交点
    % 设置起始点到下边界（y=9） (增加一个交点而不是替换)
    if ~isempty(yNewSet)
        if yNewSet(1) > 9.0
            yNewSet = [9, yNewSet];
            zNewSet = [zNewSet(1), zNewSet];
        else
            yNewSet(1) = 9;
        end
    end

    % 设置末端点到上边界（y=11） (增加一个交点而不是替换)
    if ~isempty(yNewSet)
        if yNewSet(end) < 11.0
            yNewSet = [yNewSet, 11];
            zNewSet = [zNewSet, zNewSet(end)];
        else
            yNewSet(end) = 11;
        end
    end

    % 已注释的调试代码
    % midPointUp = ceil(length(zNewSet)/2);
    % nPoint = length(zNewSet);

    %% 默认分支：其他特殊情况
    % 当所有条件都不满足时，使用备用几何约束
else
    % 调整起始点到备用上圆角
    [yNewSet(1), zNewSet(1)] = getPointOnCurveFunc([yNewSet(2), zNewSet(2)], [yNewSet(1), zNewSet(1)], upCenter2, r, 'right');

    % 调整末端点到备用上圆角
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end-1), zNewSet(end-1)], [yNewSet(end), zNewSet(end)], upCenter2, r, 'right');

end

end



