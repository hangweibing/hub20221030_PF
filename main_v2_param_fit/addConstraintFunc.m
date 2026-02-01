%% ===================================================================
%% 函数名称：addConstraintFunc
%% 功能描述：根据工件几何约束调整裂纹轮廓形状（旧版本，已被addConstraintNewSatgeFunc替代）
%% ===================================================================
function [yNewSet, zNewSet] = addConstraintFunc(yNewSet, zNewSet)

% =========================================================================
% 输入参数：
%   yNewSet:         预测的裂纹y坐标集 (21个节点)
%   zNewSet:         预测的裂纹z坐标集 (21个节点)
%
% 输出参数：
%   yNewSet:         调整后的裂纹y坐标集
%   zNewSet:         调整后的裂纹z坐标集
%
% 功能说明：
%   根据工件几何边界约束调整裂纹轮廓，确保裂纹形状符合物理约束
%   该函数已被addConstraintNewSatgeFunc替代，但保留作为历史版本
% =========================================================================

%% 初始化参数
% 查找y坐标小于7的位置（接近左边界）
minYloca = find(yNewSet < 7);

% 工件几何约束参数定义
% 下圆角几何参数 (靠近孔的下方)
downCenter1 = [6, 37.82842712];      % 下圆角圆心坐标
upCenter1 = [14, 37.82842712];       % 上圆角圆心坐标

% 其他几何约束点 (备用)
downCenter2 = [1, 58];               % 备用下圆角
upCenter2 = [19, 58];                % 备用上圆角

% 圆角半径 (mm)
r = 3;

% 已注释的调试代码
% [~, index] = sort(yNewSet);
% if(index ~= [1:11])
%     index = [1:11];
% end
% [minY, minYloca] = min(yNewSet);
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

    % 同理，处理分支2中可能过线的中间点
    over35_indices = find(zNewSet > 35);
    if ~isempty(over35_indices)
        lastZidx = over35_indices(end);
        if lastZidx > 1
            [yEnd, zEnd] = getPointOnCurveFunc([yNewSet(lastZidx), zNewSet(lastZidx)], ...
                [yNewSet(lastZidx-1), zNewSet(lastZidx-1)], ...
                upCenter1, r, 'left');
            zNewSet = [zNewSet(1:lastZidx-1), zEnd];
            yNewSet = [yNewSet(1:lastZidx-1), yEnd];
        end
    end
    %% 分支3：第二阶段 - 从下边界到圆角区域
    % 条件：无最小y位置，末端z坐标在(35, 37.83)之间，起始z坐标<30.5
    % 物理意义：裂纹扩展到圆角过渡区域
elseif(isempty(minYloca) && (35 < zNewSet(end) && zNewSet(end) < 37.82842712) && zNewSet(1) < (30 + 1e-2))
    % 调整起始点到下边界
    yNewSet(1) = interp1([zNewSet(1), zNewSet(2)], [yNewSet(1), yNewSet(2)], 30, 'linear', 'extrap'); % 基于z坐标插值调整y值到下边界
    zNewSet(1) = 30;

    % 处理上边界，可能与圆角相交的情况
    % 计算裂纹末端与上圆角的交点
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap'); % 基于y坐标插值调整z值到上边界
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    %% 分支4：第三阶段 - 从左边界到圆角区域
    % 条件：存在最小y位置，末端z坐标在(35, 37.83]之间
    % 物理意义：裂纹从左边界扩展到圆角区域
elseif((~isempty(minYloca)) && (zNewSet(end) > 35) && (zNewSet(end) <= 37.82842712))
    lastIndex = minYloca(end);
    % 计算新的起始z坐标（在y=7处）
    % newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7, 'linear', 'extrap');
    newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7); % 确保数据类型一致

    % 重构z坐标数组
    zNewSet = [newZstart, zNewSet(lastIndex+1:end)];

    % 处理上边界与圆角的交点
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    % 重构y坐标数组
    yNewSet = [7, yNewSet(lastIndex+1:end)];
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

    % 处理前缘与下圆角的相交情况，确保物理合理性
    newLeftPoint = 0;
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end
    if(newLeftPoint > 0)
        % 计算与下圆角的交点
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)] ...
            , [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        % 重构坐标数组
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    %% 分支5：第四阶段 - 上边界已过圆角，右边直线到圆角
    % 条件：存在最小y位置，末端z坐标>37.83，起始z坐标≤35
    % 物理意义：裂纹上边界已过圆角区域，右边从直线过渡到圆角
elseif((~isempty(minYloca)) && (zNewSet(end) > 37.82842712) && (zNewSet(1) <= 35))
    lastIndex = minYloca(end);
    newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7);
    zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
    yNewSet = [7, yNewSet(lastIndex+1:end)];

    % 上边界可能与圆角相交，此时需要去除圆角内部的部分
    newLeftPoint = 0;
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end
    if(newLeftPoint > 0)
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)] ...
            , [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 处理上边界
    % [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)]);
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    %% 分支6：第五阶段 - 完全在圆角区域内部
    % 条件：起始和末端z坐标都在(35, 37.83]之间
    % 物理意义：裂纹完全处于圆角区域内部
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712 && (zNewSet(1) > 35) && (zNewSet(1) <= 37.82842712))
    % 处理前缘与下圆角的相交
    newLeftPoint = 0;
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end
    if(newLeftPoint > 0)
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)] ...
            , [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 处理末端与上圆角的相交
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)], upCenter1, r, 'left');

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap'); % 基于y坐标插值调整z值到上边界
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    %% 分支7：第六阶段 - 上边界过圆角，下边界在圆角内
    % 条件：末端z坐标>37.83，起始z坐标在(35, 37.83]之间
    % 物理意义：裂纹上边界已过圆角，下边界仍在圆角内
elseif(zNewSet(end) > 37.82842712 && (zNewSet(1) > 35) && (zNewSet(1) <= 37.82842712))
    % 处理前缘与下圆角的相交
    newLeftPoint = 0;
    for i = 1:length(yNewSet)
        judge(i) = isPointInCircleFunc([yNewSet(i), zNewSet(i)], downCenter1, r);
        if(i > 1 && judge(i) == false && judge(i-1) == true)
            newLeftPoint = i - 1;
        end
    end
    if(newLeftPoint > 0)
        [yNewStart, zNewStart] = getPointOnCurveFunc([yNewSet(newLeftPoint+1), zNewSet(newLeftPoint+1)] ...
            , [yNewSet(newLeftPoint), zNewSet(newLeftPoint)], downCenter1, r, 'right');
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end
    % 已注释的复杂边界处理逻辑（备用）
    % zStartInterp = interp1([yNewSet(2), yNewSet(3)], [zNewSet(2), zNewSet(3)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % if(zNewSet(1) > zStartInterp)
    %     zNewSet(1) = zStartInterp;
    % else
    %     zNewSet(1) = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % end

    % 处理上边界
    zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap'); % 基于y坐标插值调整z值到上边界
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');
    %% 分支8：第七阶段 - 完全在上边界区域
    % 条件：起始z坐标在(37.83, 50)之间，末端z坐标<50
    % 物理意义：裂纹完全处于上边界直线区域
elseif(zNewSet(1) > 37.82842712 && (zNewSet(end) < 50) && (zNewSet(1) > 37.82842712) && (zNewSet(1) < 50))
    % 处理下边界
    zNewSet(1) = interp1([yNewSet(2), yNewSet(3)], [zNewSet(2), zNewSet(3)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');

    % 已注释的复杂处理逻辑（备用）
    % midPointDown = floor(length(zNewSet)/2);
    % nPoint = length(zNewSet);
    % [minDown, minDownLoca] = min(zNewSet(1:midPointDown));
    % ... 复杂的边界调整逻辑 ...

    % 处理上边界
    yNewSet(1) = getEdgeYbyZFunc(zNewSet(1), 'down');
    zNewSet(end) = interp1([yNewSet(end-2), yNewSet(end-1)], [zNewSet(end-2), zNewSet(end-1)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap');
    midPointUp = ceil(length(zNewSet)/2);

    % 已注释的复杂处理逻辑（备用）
    % [minUp, minUpLoca] = min(zNewSet(midPointUp:end));
    % ... 复杂的边界调整逻辑 ...

    zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap'); % 基于y坐标插值调整z值到上边界
    yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

    %% 分支9：其他特殊情况
    % 条件：不满足上述任何条件的情况
    % 物理意义：使用备用几何约束处理
    % 使用备用圆角处理（upCenter2）
    [yNewSet(1), zNewSet(1)] = getPointOnCurveFunc([yNewSet(2), zNewSet(2)], [yNewSet(1), zNewSet(1)], upCenter2, r, 'right');
    [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end-1), zNewSet(end-1)], [yNewSet(end), zNewSet(end)], upCenter2, r, 'right');
end

%% ===================================================================
%% 辅助函数定义
%% ===================================================================

%% ===================================================================
%% 函数名称：getPointOnCurveFunc
%% 功能描述：计算直线与圆的交点
%% ===================================================================
    function [x, y] = getPointOnCurveFunc(point1, point2, center, r, location)
    end

%% GETPOINTONCURVESCRIPT 计算直线与圆的交点
%
% 输入参数：
%   point1, point2: 直线上的两个点 [x, y]
%   center: 圆心坐标 [a, b]
%   r: 圆半径
%   location: 选择交点位置 ('left' 或 'right')
%
% 输出参数：
%   x, y: 选定的交点坐标
%
% 数学原理：
% 圆方程: (x-a)² + (y-b)² = r²
% 直线方程: c*x + d*y = e

% 提取圆心坐标
a = center(1);  % 圆心x坐标
b = center(2);  % 圆心y坐标
% r = 3;  % 圆半径（已通过参数传入）

% 提取直线上的两个点
x1 = point1(1);
y1 = point1(2);
x2 = point2(1);
y2 = point2(2);

% 计算直线方程系数
% 直线方程: c*x + d*y = e
c = y2 - y1;     % x系数
d = x1 - x2;     % y系数
e = x1*y2 - x2*y1;  % 常数项
% ===================================================================
% 求解直线与圆的交点（联立方程组）
% 圆方程: (x-a)² + (y-b)² = r²
% 直线方程: c*x + d*y = e
%
% 使用代数方法求解，得到两个交点 (xc1,yc1) 和 (xc2,yc2)
% ===================================================================

% 第一个交点坐标（解方程组得到）
xc1 = (1/2) .* (c.^2 + d.^2).^(-1) .* ((-2).*b.*c.*d + 2.*a.*d.^2 + 2.*c.*e + ...
    (-1).*((2.*b.*c.*d + (-2).*a.*d.^2 + (-2).*c.*e).^2 + (-4).*(c.^2 + ...
    d.^2).*(a.^2.*d.^2 + b.^2.*d.^2 + (-2).*b.*d.*e + e.^2 + (-1).* ...
    d.^2.*r.^2)).^(1/2));

yc1 = -(d.^(-1) .* ((-1).*b.*c.^2.*d.*(c.^2 + d.^2).^(-1) + a.*c.*d.^2.*( ...
    c.^2 + d.^2).^(-1) + (-1).*e + c.^2.*(c.^2 + d.^2).^(-1).*e + (-1/2).* ...
    c.*(c.^2 + d.^2).^(-1).*((2.*b.*c.*d + (-2).*a.*d.^2 + (-2).*c.*e) ...
    .^2 + (-4).*(c.^2 + d.^2).*(a.^2.*d.^2 + b.^2.*d.^2 + (-2).*b.*d.*e + ...
    e.^2 + (-1).*d.^2.*r.^2)).^(1/2)));

% 第二个交点坐标（解方程组得到）
xc2 = (1/2).*(c.^2 + d.^2).^(-1).*((-2).*b.*c.*d + 2.*a.*d.^2 + 2.*c.*e + ...
    ((2.*b.*c.*d + (-2).*a.*d.^2 + (-2).*c.*e).^2 + (-4).*(c.^2 + d.^2) ...
    .*(a.^2.*d.^2 + b.^2.*d.^2 + (-2).*b.*d.*e + e.^2 + (-1).*d.^2.* ...
    r.^2)).^(1/2));

yc2 = (-1).*d.^(-1).*((-1).*b.*c.^2.*d.*(c.^2 + d.^2).^(-1) + a.*c.* ...
    d.^2.*(c.^2 + d.^2).^(-1) + (-1).*e + c.^2.*(c.^2 + d.^2).^(-1).*e + ( ...
    1/2).*c.*(c.^2 + d.^2).^(-1).*((2.*b.*c.*d + (-2).*a.*d.^2 + (-2) ...
    .*c.*e).^2 + (-4).*(c.^2 + d.^2).*(a.^2.*d.^2 + b.^2.*d.^2 + (-2).* ...
    b.*d.*e + e.^2 + (-1).*d.^2.*r.^2)).^(1/2));
% ===================================================================
% 根据location参数选择合适的交点
% ===================================================================

if strcmp(location, 'left')
    % 选择x坐标较小的交点（左侧）
    if(xc1 < xc2)
        x = xc1;
        y = yc1;
    else
        x = xc2;
        y = yc2;
    end
elseif strcmp(location, 'right')
    % 选择x坐标较大的交点（右侧）
    if(xc1 < xc2)
        x = xc2;
        y = yc2;
    else
        x = xc1;
        y = yc1;
    end
else
    % 参数错误提示
    error('位置参数错误，请使用''left''或''right''')
end

%% ===================================================================
%% 函数名称：isPointInCircleFunc
%% 功能描述：判断点是否在圆内
%% ===================================================================
    function judge = isPointInCircleFunc(point, center, r)
    end

%% ISPOINTINCIRCLEFUNC 判断点是否在圆内
%
% 输入参数：
%   point: 待判断点坐标 [x, y]
%   center: 圆心坐标 [a, b]
%   r: 圆半径
%
% 输出参数：
%   judge: 布尔值，true表示点在圆内（含边界），false表示在圆外

% 计算点到圆心的距离
distance = sqrt((point(1) - center(1)).^2 + (point(2) - center(2)).^2);

% 判断距离是否小于等于半径
if(distance <= r)
    judge = true;   % 点在圆内（含边界）
else
    judge = false;  % 点在圆外
end
end

