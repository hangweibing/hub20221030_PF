%% ===================================================================
%% 函数名称：addConstraintNewSatgeFunc
%% 功能描述：根据工件几何约束调整裂纹轮廓形状
%% ===================================================================
function [yNewSet, zNewSet, splitted] = addConstraintNewSatgeFunc(yNewSet, zNewSet, logCstar, gamma)

% =========================================================================
% 输入参数：
%   yNewSet:         预测的裂纹y坐标集 (21个节点)
%   zNewSet:         预测的裂纹z坐标集 (21个节点)
%   logCstar:        Paris定律参数 logC* (用于调试可视化)
%   gamma:           Paris定律参数 γ (用于调试可视化)
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

%% 保存原始输入数据（用于调试对比）
yNewSet_original = yNewSet;          % 处理前的y坐标（完整21节点）
zNewSet_original = zNewSet;          % 处理前的z坐标（完整21节点）

%% 数据有效性检查阈值（统一参数控制）
MIN_CRACK_LENGTH = 1;                % 最小合理裂纹长度（mm）
MAX_CRACK_LENGTH = 23;               % 最大合理裂纹长度（mm）

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

%% 分支1：第一阶段 - 直线到上边界
% 条件：无最小y位置且末端z坐标≤35
% 物理意义：裂纹处于初始阶段，主要在直线区域扩展
if(isempty(minYloca) && zNewSet(end) <= 35)
    % 首先调整起始点到下边界
    yNewSet(1) = interp1([zNewSet(1), zNewSet(2)], [yNewSet(1), yNewSet(2)], 30, 'linear', 'extrap'); % 基于z坐标插值调整y值到下边界
    zNewSet(1) = 30;  % 固定起始z坐标到下边界

    % 调整末端点到上边界
    zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], 13, 'linear', 'extrap'); % 基于y坐标插值调整z值
    yNewSet(end) = 13;  % 固定末端y坐标到上边界

    % 保存几何调整后的状态（用于可视化）
    yAfterClean = yNewSet;
    zAfterClean = zNewSet;

    % 检查点：分支1结束
    checkCrackLength(yNewSet, zNewSet, '分支1-直线到上边界', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma, yAfterClean, zAfterClean);

%% 分支2：第一阶段实际分支 - 包含内部拐点，直线到直线
% 条件：存在最小y位置且末端z坐标≤35
% 物理意义：裂纹有内部拐点，从下边界开始向上扩展
elseif(~isempty(minYloca) && zNewSet(end) <= 35)
    lastIndex = minYloca(end);  % 最后一个最小y位置的索引
    % 计算新的起始z坐标（在y=7处）
    newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7, 'linear', 'extrap');

    % 重构z坐标数组，从新的起始点开始
    zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
    % 调整末端点到上边界
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    % 重构y坐标数组，从y=7开始
    yNewSet = [7, yNewSet(lastIndex+1:end)];

    % 保存几何调整后的状态（用于可视化）
    yAfterClean = yNewSet;
    zAfterClean = zNewSet;

    % 检查点：分支2结束
    checkCrackLength(yNewSet, zNewSet, '分支2-内部拐点', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma, yAfterClean, zAfterClean);

%% 分支3：第二阶段 - 从下边界到圆角区域
% 条件：无最小y位置，末端z坐标在(35, 37.83)之间，起始z坐标<30.5
% 物理意义：裂纹扩展到圆角过渡区域
elseif(isempty(minYloca) && (35 < zNewSet(end) && zNewSet(end) < 37.82842712) && zNewSet(1) < 30.5)
    % 调整起始点到下边界
    yNewSet(1) = interp1([zNewSet(1), zNewSet(2)], [yNewSet(1), yNewSet(2)], 30, 'linear', 'extrap');
    zNewSet(1) = 30;

    % 处理上边界，可能与圆角相交的情况
    % 计算裂纹末端与上圆角的交点
    [y_tmp, z_tmp, COMPLEX_SOULTION] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    if COMPLEX_SOULTION  % 如果没有有效交点
        yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');  % 使用上边界函数
    else
        yNewSet(end) = y_tmp;  % 使用计算的交点
        zNewSet(end) = z_tmp;
    end

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap');
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

    % 保存几何调整后的状态（用于可视化）
    yAfterClean = yNewSet;
    zAfterClean = zNewSet;

    % 检查点：分支3结束
    checkCrackLength(yNewSet, zNewSet, '分支3-圆角过渡', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma, yAfterClean, zAfterClean);

%% 分支4：第三阶段 - 从圆角到上边界（判断交点）
% 条件：(存在最小y位置 或 起始z坐标在(30,35)之间) 且 末端z坐标在(35,37.83]之间
% 物理意义：裂纹在圆角区域内扩展，可能与下圆角相交
elseif (~isempty(minYloca) || ((zNewSet(1)>30) && (zNewSet(1)<35))) && (zNewSet(end)>35) && (zNewSet(end)<=37.82842712)
    % 处理起始点调整
    if ~isempty(minYloca)
        lastIndex = minYloca(end);
        % 计算在y=7处的z坐标
        newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7);
        % 重构坐标数组
        zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
        yNewSet = [7, yNewSet(lastIndex+1:end)];
    else
        % 计算新的起始z坐标
        newZstart = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], 7, 'linear', 'extrap');
        zNewSet = [newZstart, zNewSet(2:end)];
        yNewSet = [7, yNewSet(2:end)];
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    % 检查是否与下圆角相交，如果是则进行前缘分裂处理
    [yNewSet, zNewSet, splitted] = resplit_front(yNewSet, zNewSet, downCenter1, r);
    
    % 检查点：分支4结束
    checkCrackLength(yNewSet, zNewSet, '分支4-圆角区域', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma);
    
%% 分支5：第四阶段 - 圆角到圆角
% 条件：起始和末端z坐标都在(35,37.83]之间
% 物理意义：裂纹完全在圆角区域内扩展
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712 && zNewSet(1) > 35 && zNewSet(1) <= 37.82842712)
    newLeftPoint = 0;  % 新的左边界点索引

    % 检查哪些点在下圆角内
    judge = zeros(1, n_reg_point);
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        % 寻找从圆内到圆外的过渡点
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end

    % 如果找到过渡点，调整起始点到圆角边界
    if(newLeftPoint > 0)
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)], ...
            [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap');
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    
    % 检查点：分支5结束
    checkCrackLength(yNewSet, zNewSet, '分支5-圆角到圆角', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma);

%% 分支6：第五阶段 - 从圆角到直线区域（判断交点）
% 条件：末端z坐标在(35,37.83]之间，起始z坐标≥37.83
% 物理意义：裂纹从圆角区域扩展到直线区域
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712) && (zNewSet(1) >= 37.82842712)
    % 首先清理多余的直线部分
    min_y_loca = find(yNewSet < 9);  % 查找y坐标小于9的位置（接近左边界）
    if ~isempty(min_y_loca)
        % 保留从最后一个小y值位置到末尾的部分
        zNewSet = zNewSet(min_y_loca(end):end);
        yNewSet = yNewSet(min_y_loca(end):end);
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');
    
    % 检查点：分支6结束
    checkCrackLength(yNewSet, zNewSet, '分支6-圆角到直线', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma);
    
%% 分支7：第六阶段 - 从直线到上边界（判断交点）
% 条件：末端z坐标>37.83，起始z坐标≤35
% 物理意义：裂纹从下部直线区域扩展到上边界，可能与下圆角相交
elseif(zNewSet(end) > 37.82842712) && (zNewSet(1) <= 35)
    % 处理起始点调整到左边界
    if ~isempty(minYloca)
        lastIndex = minYloca(end);
        newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7);
        zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
        yNewSet = [7, yNewSet(lastIndex+1:end)];
    else
        newZstart = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], 7, 'linear', 'extrap');
        zNewSet = [newZstart, zNewSet(2:end)];
        yNewSet = [7, yNewSet(2:end)];
    end

    % 调整末端点到上边界
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

    % 检查是否与下圆角相交
    [yNewSet, zNewSet, splitted] = resplit_front(yNewSet, zNewSet, downCenter1, r);

    % 已注释的备用处理方法
    % [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)]);
    
    % 检查点：分支7结束
    checkCrackLength(yNewSet, zNewSet, '分支7-直线到上边界', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma);

%% 分支8：第七阶段 - 从直线到圆角
% 条件：末端z坐标>37.83，起始z坐标在(35,37.83]之间
% 物理意义：裂纹从直线区域进入圆角区域
elseif(zNewSet(end) > 37.82842712 && zNewSet(1) > 35 && zNewSet(1) <= 37.82842712)
    newLeftPoint = 0;  % 新的左边界点

    % 检查点与下圆角的相对位置
    judge = zeros(1, n_reg_point);
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        % 寻找从圆内到圆外的过渡点
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end


    % 调整起始点到圆角边界
    if(newLeftPoint > 0)
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)], ...
            [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 保存几何调整后的状态（用于可视化）
    yAfterClean = yNewSet;
    zAfterClean = zNewSet;

    % 调整末端点到上边界（使用插值方法）
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');



    % 已注释的备用起始点处理方法
    % zStartInterp = interp1([yNewSet(2), yNewSet(3)], [zNewSet(2), zNewSet(3)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % if(zNewSet(1) > zStartInterp)
    %     zNewSet(1) = zStartInterp;
    % else
    %     zNewSet(1) = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % end

    % 检查点：分支8结束
    checkCrackLength(yNewSet, zNewSet, '分支8-直线到圆角', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma, yAfterClean, zAfterClean);
    
%% 分支9：第八阶段 - 上边界区域（直线到直线）
% 条件：起始和末端z坐标都>37.83
% 物理意义：裂纹在上边界直线区域内扩展
elseif zNewSet(end) > 37.82842712 && zNewSet(1) > 37.82842712  % 20220915 修改：上边界区域的处理
    % 清理下边界多余部分
    min_y_loca = find(yNewSet < 9);  % 查找在下边界(y<9)外的点
    if ~isempty(min_y_loca)
        % 保留从最后一个小y值位置到末尾的部分
        zNewSet = zNewSet(min_y_loca(end):end);
        yNewSet = yNewSet(min_y_loca(end):end);
    end

    % 清理上边界多余部分 (220915新增)
    max_y_loca = find(yNewSet > 11);  % 查找超过上边界(y>11)的点
    if ~isempty(max_y_loca)
        % 保留从开始到第一个大y值位置的部分
        zNewSet = zNewSet(1:max_y_loca(end));
        yNewSet = yNewSet(1:max_y_loca(end));
    end


    % 调整起始点到下边界
    zNewSet(1) = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], getEdgeYbyZFunc(zNewSet(1), 'down'), 'linear', 'extrap');
    yNewSet(1) = getEdgeYbyZFunc(zNewSet(1), 'down');

    
    % 保存边界清理后的状态（用于可视化）
    yAfterClean = yNewSet;
    zAfterClean = zNewSet;

    % 调整末端点到上边界（直接投影方法）
    % 说明：当末端几个点的y坐标几乎不变时，插值方法会产生极大偏差
    % 因此采用简单的垂直投影方法：保持z坐标不变，只将y坐标投影到上边界
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    % zNewSet(end) = zNewSet(end);  % z坐标保持不变（可省略）
    
    % 已注释的插值方法（当末端y坐标几乎不变时会产生偏差）
    % zNewSet(end) = interp1([yNewSet(end-2), yNewSet(end-1)], [zNewSet(end-2), zNewSet(end-1)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

    % 已注释的调试代码
    % midPointUp = ceil(length(zNewSet)/2);
    % nPoint = length(zNewSet);
    
    % 检查点：分支9结束
    checkCrackLength(yNewSet, zNewSet, '分支9-上边界区域', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma, yAfterClean, zAfterClean);

%% 默认分支：其他特殊情况
% 当所有条件都不满足时，使用备用几何约束
else
    % 调整起始点到备用上圆角
    [yNewSet(1), zNewSet(1)] = getPointOnCurveFunc([yNewSet(2), zNewSet(2)], [yNewSet(1), zNewSet(1)], upCenter2, r, 'right');

    % 调整末端点到备用上圆角
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end-1), zNewSet(end-1)], [yNewSet(end), zNewSet(end)], upCenter2, r, 'right');
    
    % 检查点：默认分支结束
    checkCrackLength(yNewSet, zNewSet, '默认分支-备用约束', MIN_CRACK_LENGTH, MAX_CRACK_LENGTH, yNewSet_original, zNewSet_original, logCstar, gamma);
end

end

%% ===================================================================
%% 内部函数：数据有效性检查
%% ===================================================================
function checkCrackLength(yNewSet, zNewSet, branch_name, min_len, max_len, yOriginal, zOriginal, logCstar, gamma, yAfterClean, zAfterClean)
    % 检查裂纹长度是否在合理范围内
    % 输入参数：
    %   yNewSet, zNewSet: 处理后的坐标（可能节点数<21）
    %   yOriginal, zOriginal: 处理前的原始坐标（完整21节点）
    %   logCstar, gamma: Paris定律参数（用于可视化标注）
    %   yAfterClean, zAfterClean: 边界清理后的坐标（可选，仅分支9使用）
    
    % 判断是否提供了边界清理后数据
    hasCleanData = (nargin >= 11) && ~isempty(yAfterClean) && ~isempty(zAfterClean);
    
    upcrack_z = zNewSet(end);
    crack_length = upcrack_z - 30;
    n_points_before = length(yOriginal);  % 处理前节点数（通常是21）
    n_points_after = length(yNewSet);     % 处理后节点数（可能<21）
    if hasCleanData
        n_points_clean = length(yAfterClean);  % 边界清理后节点数
    end
    
    % 自动计算显示范围（增加1mm边距）
    if hasCleanData
        all_z = [real(zOriginal(:)); real(zAfterClean(:)); real(zNewSet(:))];
    else
        all_z = [real(zOriginal(:)); real(zNewSet(:))];
    end
    z_min = min(all_z) - 1;
    z_max = max(all_z) + 1;
    
    % 检查是否为复数
    if ~isreal(upcrack_z) || ~isreal(yNewSet) || ~isreal(zNewSet)
        % 可视化异常坐标并保存为文件（parfor兼容）- 三步对比图
        if hasCleanData
            fig = figure('Visible', 'off', 'Position', [100, 100, 1800, 400]);
        else
            fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);
        end
        
        % 左图：处理前（原始输入，21节点）
        if hasCleanData
            subplot(1, 3, 1);
        else
            subplot(1, 2, 1);
        end
        plot(zOriginal, yOriginal, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        plot(zOriginal(end), yOriginal(end), 'g^', 'MarkerSize', 8, 'LineWidth', 1.5);

        % 绘制物理边界
        z_boundary = 30:0.1:50;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([30, 50]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        title(sprintf('处理前（节点数=%d）', n_points_before));
        grid on;
        legend('原始节点', '末端点', 'Location', 'best');

        % 中图：边界清理后（仅当有清理数据时显示）
        if hasCleanData
            subplot(1, 3, 2);
            plot(zAfterClean, yAfterClean, 'mo-', 'LineWidth', 1.5, 'MarkerSize', 4);
            hold on;
            plot(zAfterClean(end), yAfterClean(end), 'c^', 'MarkerSize', 8, 'LineWidth', 1.5);

            % 绘制物理边界
            z_boundary = 30:0.1:50;
            y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
            y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
            plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
            plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

            xlim([30, 50]);
            ylim([6, 15]);
            xlabel('z坐标 (mm)');
            ylabel('y坐标 (mm)');
            title(sprintf('边界清理后（节点数=%d）', n_points_clean));
            grid on;
            legend('清理后节点', '末端点', 'Location', 'best');
        end

        % 右图：最终处理后（出现复数，节点数可能变化）
        if hasCleanData
            subplot(1, 3, 3);
        else
            subplot(1, 2, 2);
        end
        plot(real(zNewSet), real(yNewSet), 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4);
        hold on;
        plot(real(zNewSet(end)), real(yNewSet(end)), 'mx', 'MarkerSize', 10, 'LineWidth', 2);
        
        % 绘制物理边界
        z_boundary = z_min:0.1:z_max;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([z_min, z_max]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        title(sprintf('处理后-异常复数值（节点数=%d）', n_points_after));
        grid on;
        legend('处理后节点', '末端点（异常）', 'Location', 'best');
        
        sgtitle(sprintf('异常检测：复数值 - %s\n粒子参数: logC*=%.3f, γ=%.3f', branch_name, logCstar, gamma), 'FontSize', 14, 'FontWeight', 'bold');
        
        % 保存图形到结果文件夹 (使用 saveas 以确保更好的 SVG 兼容性)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('结果/异常检测_复数值_%s_%s.svg', strrep(branch_name, '-', '_'), timestamp);
        saveas(fig, filename); 
        close(fig);
        
        error(['[addConstraintNewSatgeFunc-%s] 检测到复数值！\n' ...
               '  节点数：处理前%d → 处理后%d\n' ...
               '  上表面z坐标: %.6f\n' ...
               '  可视化已保存(SVG矢量图): %s\n' ...
               '  请检查几何插值或边界处理。'], branch_name, n_points_before, n_points_after, upcrack_z, filename);
    end
    
    % 检查是否为异常小值
    if crack_length < min_len
        % 可视化异常坐标并保存为文件（parfor兼容）- 三步对比图
        if hasCleanData
            fig = figure('Visible', 'off', 'Position', [100, 100, 1800, 400]);
        else
            fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);
        end

        % 左图：处理前（原始输入，21节点）
        if hasCleanData
            subplot(1, 3, 1);
        else
            subplot(1, 2, 1);
        end
        plot(zOriginal, yOriginal, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        plot(zOriginal(end), yOriginal(end), 'g^', 'MarkerSize', 8, 'LineWidth', 1.5);

        % 绘制物理边界
        z_boundary = 30:0.1:50;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([30, 50]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        crack_len_before = zOriginal(end) - 30;
        title(sprintf('处理前（节点数=%d）\n裂纹长度=%.2f mm', n_points_before, crack_len_before));
        grid on;
        legend('原始节点', '末端点', 'Location', 'best');

        % 中图：边界清理后（仅当有清理数据时显示）
        if hasCleanData
            subplot(1, 3, 2);
            plot(zAfterClean, yAfterClean, 'mo-', 'LineWidth', 1.5, 'MarkerSize', 4);
            hold on;
            plot(zAfterClean(end), yAfterClean(end), 'c^', 'MarkerSize', 8, 'LineWidth', 1.5);

            % 绘制物理边界
            z_boundary = 30:0.1:50;
            y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
            y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
            plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
            plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

            xlim([30, 50]);
            ylim([6, 15]);
            xlabel('z坐标 (mm)');
            ylabel('y坐标 (mm)');
            crack_len_clean = zAfterClean(end) - 30;
            title(sprintf('边界清理后（节点数=%d）\n裂纹长度=%.2f mm', n_points_clean, crack_len_clean));
            grid on;
            legend('清理后节点', '末端点', 'Location', 'best');
        end

        % 右图：最终处理后（异常小值）
        if hasCleanData
            subplot(1, 3, 3);
        else
            subplot(1, 2, 2);
        end
        plot(zNewSet, yNewSet, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4);
        hold on;
        plot(zNewSet(end), yNewSet(end), 'mx', 'MarkerSize', 10, 'LineWidth', 2);
        
        % 绘制物理边界
        z_boundary = z_min:0.1:z_max;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([z_min, z_max]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        title(sprintf('处理后-异常小值（节点数=%d）\n裂纹长度=%.2f mm < %.1f mm', n_points_after, crack_length, min_len));
        grid on;
        legend('处理后节点', '末端点（异常）', 'Location', 'best');
        
        sgtitle(sprintf('异常检测：裂纹长度过小 - %s\n粒子参数: logC*=%.3f, γ=%.3f', branch_name, logCstar, gamma), 'FontSize', 14, 'FontWeight', 'bold');
        
        % 保存图形到结果文件夹 (使用 saveas 以确保更好的 SVG 兼容性)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('结果/异常检测_异常小值_%s_%s.svg', strrep(branch_name, '-', '_'), timestamp);
        saveas(fig, filename);
        close(fig);
        
        error(['[addConstraintNewSatgeFunc-%s] 检测到异常小值！\n' ...
               '  节点数：处理前%d → 处理后%d\n' ...
               '  上表面z坐标: %.3f mm\n' ...
               '  裂纹长度: 处理前%.3f mm → 处理后%.3f mm (小于阈值 %.1f mm)\n' ...
               '  可视化已保存(SVG矢量图): %s\n' ...
               '  zNewSet范围: [%.3f, %.3f]\n' ...
               '  yNewSet范围: [%.3f, %.3f]'], ...
               branch_name, n_points_before, n_points_after, upcrack_z, ...
               crack_len_before, crack_length, min_len, filename, ...
               min(zNewSet), max(zNewSet), min(yNewSet), max(yNewSet));
    end
    
    % 检查是否超过合理阈值
    if crack_length > max_len
        % 可视化异常坐标并保存为文件（parfor兼容）- 三步对比图
        if hasCleanData
            fig = figure('Visible', 'off', 'Position', [100, 100, 1800, 400]);
        else
            fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);
        end

        % 左图：处理前（原始输入，21节点）
        if hasCleanData
            subplot(1, 3, 1);
        else
            subplot(1, 2, 1);
        end
        plot(zOriginal, yOriginal, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 3);
        hold on;
        plot(zOriginal(end), yOriginal(end), 'g^', 'MarkerSize', 8, 'LineWidth', 1.5);

        % 绘制物理边界
        z_boundary = 30:0.1:50;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([30, 50]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        crack_len_before = zOriginal(end) - 30;
        title(sprintf('处理前（节点数=%d）\n裂纹长度=%.2f mm', n_points_before, crack_len_before));
        grid on;
        legend('原始节点', '末端点', 'Location', 'best');

        % 中图：边界清理后（仅当有清理数据时显示）
        if hasCleanData
            subplot(1, 3, 2);
            plot(zAfterClean, yAfterClean, 'mo-', 'LineWidth', 1.5, 'MarkerSize', 4);
            hold on;
            plot(zAfterClean(end), yAfterClean(end), 'c^', 'MarkerSize', 8, 'LineWidth', 1.5);

            % 绘制物理边界
            z_boundary = 30:0.1:50;
            y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
            y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
            plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
            plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

            xlim([30, 50]);
            ylim([6, 15]);
            xlabel('z坐标 (mm)');
            ylabel('y坐标 (mm)');
            crack_len_clean = zAfterClean(end) - 30;
            title(sprintf('边界清理后（节点数=%d）\n裂纹长度=%.2f mm', n_points_clean, crack_len_clean));
            grid on;
            legend('清理后节点', '末端点', 'Location', 'best');
        end

        % 右图：最终处理后（异常大值）
        if hasCleanData
            subplot(1, 3, 3);
        else
            subplot(1, 2, 2);
        end
        plot(zNewSet, yNewSet, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 4);
        hold on;
        plot(zNewSet(end), yNewSet(end), 'mx', 'MarkerSize', 10, 'LineWidth', 2);
        
        % 绘制物理边界
        z_boundary = z_min:0.1:z_max;
        y_upper = arrayfun(@(z) getEdgeYbyZFunc(z, 'up'), z_boundary);
        y_lower = arrayfun(@(z) getEdgeYbyZFunc(z, 'down'), z_boundary);
        plot(z_boundary, y_upper, 'k-', 'LineWidth', 0.8);
        plot(z_boundary, y_lower, 'k-', 'LineWidth', 0.8);

        xlim([z_min, z_max]);
        ylim([6, 15]);
        xlabel('z坐标 (mm)');
        ylabel('y坐标 (mm)');
        title(sprintf('处理后-异常大值（节点数=%d）\n裂纹长度=%.2f mm > %.1f mm', n_points_after, crack_length, max_len));
        grid on;
        legend('处理后节点', '末端点（异常）', 'Location', 'best');
        
        sgtitle(sprintf('异常检测：裂纹长度过大 - %s\n粒子参数: logC*=%.3f, γ=%.3f', branch_name, logCstar, gamma), 'FontSize', 14, 'FontWeight', 'bold');
        
        % 保存图形到结果文件夹 (使用 saveas 以确保更好的 SVG 兼容性)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('结果/异常检测_异常大值_%s_%s.svg', strrep(branch_name, '-', '_'), timestamp);
        saveas(fig, filename);
        close(fig);
        
        error(['[addConstraintNewSatgeFunc-%s] 检测到异常大值！\n' ...
               '  节点数：处理前%d → 处理后%d\n' ...
               '  上表面z坐标: %.3f mm\n' ...
               '  裂纹长度: 处理前%.3f mm → 处理后%.3f mm (超过阈值 %.1f mm)\n' ...
               '  可视化已保存(SVG矢量图): %s\n' ...
               '  zNewSet范围: [%.3f, %.3f]\n' ...
               '  yNewSet范围: [%.3f, %.3f]'], ...
               branch_name, n_points_before, n_points_after, upcrack_z, ...
               crack_len_before, crack_length, max_len, filename, ...
               min(zNewSet), max(zNewSet), min(yNewSet), max(yNewSet));
    end
end



