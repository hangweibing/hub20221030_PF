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
    % yNewSet(1) = interp1([zNewSet(1), zNewSet(2)], [yNewSet(1), yNewSet(2)], 30, 'linear', 'extrap'); % 基于z坐标插值调整y值到下边界
    zNewSet(1) = 30;  % 直接投影坐标，不插值

    % 调整末端点到上边界
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], 13, 'linear', 'extrap'); % 基于y坐标插值调整z值
    yNewSet(end) = 13;  % 固定末端y坐标到上边界

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 1');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

%% 分支2：
elseif(~isempty(minYloca) && zNewSet(end) <= 35)
    lastIndex = minYloca(end);  % 最后一个最小y位置的索引
    % 计算新的起始z坐标（在y=7处）
    newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7, 'linear', 'extrap');

    % 重构z坐标数组，从新的起始点开始
    zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
    % 调整末端点到上边界
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));
    % 重构y坐标数组，从y=7开始
    yNewSet = [7, yNewSet(lastIndex+1:end)];

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 2');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

%% 分支3：第二阶段
elseif(isempty(minYloca) && (35 < zNewSet(end) && zNewSet(end) < 37.82842712) && zNewSet(1) < 30.5)
    % 调整起始点到下边界
    % yNewSet(1) = interp1([zNewSet(1), zNewSet(2)], [yNewSet(1), yNewSet(2)], 30, 'linear', 'extrap');
    zNewSet(1) = 30;  % 直接投影坐标，不插值

    % 处理上边界，可能与圆角相交的情况
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 3');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap');
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up');

%% 分支4：
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
        % newZstart = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], 7, 'linear', 'extrap');
        yNewSet = [7, yNewSet(2:end)];% 直接投影，不插值
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    % 检查是否与下圆角相交，如果是则进行前缘分裂处理
    [yNewSet, zNewSet, splitted] = resplit_front(yNewSet, zNewSet, downCenter1, r);

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 4');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end
    
%% 分支5：
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
        [yNewStart, zNewStart] = getEdgeYbyZFunc(zNewSet(newLeftPoint), 'down', yNewSet(newLeftPoint));
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 5');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

    % 已注释的备用处理方法
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end),'up'), 'linear', 'extrap');
    % yNewSet(end) = getEdgeYbyZFunc(zNewSet(end), 'up'); 

%% 分支6：
elseif(zNewSet(end) > 35 && zNewSet(end) <= 37.82842712) && (zNewSet(1) >= 37.82842712)
    % 首先清理多余的直线部分
    min_y_loca = find(yNewSet < 9);  % 查找y坐标小于9的位置（接近左边界）
    if ~isempty(min_y_loca)
        % 保留从最后一个小y值位置到末尾的部分
        zNewSet = zNewSet(min_y_loca(end):end);
        yNewSet = yNewSet(min_y_loca(end):end);
    end

    % 调整末端点到上圆角
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 6');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end
    
%% 分支7：
elseif(zNewSet(end) > 37.82842712) && (zNewSet(1) <= 35)
    % 处理起始点调整到左边界
    if ~isempty(minYloca)
        lastIndex = minYloca(end);
        newZstart = interp1([yNewSet(lastIndex), yNewSet(lastIndex+1)], [zNewSet(lastIndex), zNewSet(lastIndex+1)], 7);
        zNewSet = [newZstart, zNewSet(lastIndex+1:end)];
        yNewSet = [7, yNewSet(lastIndex+1:end)];
    else
        % newZstart = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], 7, 'linear', 'extrap');
        yNewSet = [7, yNewSet(2:end)];% 直接投影，不插值
    end

    % 调整末端点到上边界
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    % 检查是否与下圆角相交
    [yNewSet, zNewSet, splitted] = resplit_front(yNewSet, zNewSet, downCenter1, r);

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 7');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

    % 已注释的备用处理方法
    % [yNewSet(end), zNewSet(end)] = getPointOnCurveFunc([yNewSet(end), zNewSet(end)], [yNewSet(end-1), zNewSet(end-1)]);
    
    

%% 分支8：
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
        [yNewStart, zNewStart] = getEdgeYbyZFunc(zNewSet(newLeftPoint), 'down', yNewSet(newLeftPoint));
        zNewSet = [zNewStart, zNewSet(newLeftPoint+1:end)];
        yNewSet = [yNewStart, yNewSet(newLeftPoint+1:end)];
    end

    % 调整末端点到上边界（使用插值方法）
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 8');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end

    % 已注释的备用起始点处理方法
    % zStartInterp = interp1([yNewSet(2), yNewSet(3)], [zNewSet(2), zNewSet(3)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % if(zNewSet(1) > zStartInterp)
    %     zNewSet(1) = zStartInterp;
    % else
    %     zNewSet(1) = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], getEdgeYbyZFunc(zNewSet(1),'down'), 'linear', 'extrap');
    % end
    
%% 分支9：
elseif zNewSet(end) > 37.82842712 && zNewSet(1) > 37.82842712  % 20220915 修改：上边界区域的处理
    % 清理左边界多余部分
    min_y_loca = find(yNewSet < 9);  % 查找接近左边界(y<9)的点
    if ~isempty(min_y_loca)
        % 保留从最后一个小y值位置到末尾的部分
        zNewSet = zNewSet(min_y_loca(end):end);
        yNewSet = yNewSet(min_y_loca(end):end);
    end

    % 清理右边界多余部分 (220915新增)
    max_y_loca = find(yNewSet > 11);  % 查找超过右边界(y>11)的点
    if ~isempty(max_y_loca)
        % 保留从开始到第一个大y值位置的部分
        zNewSet = zNewSet(1:max_y_loca(end));
        yNewSet = yNewSet(1:max_y_loca(end));
    end

    % 调整起始点到下边界
    % zNewSet(1) = interp1([yNewSet(1), yNewSet(2)], [zNewSet(1), zNewSet(2)], getEdgeYbyZFunc(zNewSet(1), 'down'), 'linear', 'extrap');
    [yNewSet(1), zNewSet(1)] = getEdgeYbyZFunc(zNewSet(1), 'down', yNewSet(1));

    % 调整末端点到上边界（两种插值方法）
    % zNewSet(end) = interp1([yNewSet(end-2), yNewSet(end-1)], [zNewSet(end-2), zNewSet(end-1)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    % 最终调整末端点
    % zNewSet(end) = interp1([yNewSet(end-1), yNewSet(end)], [zNewSet(end-1), zNewSet(end)], getEdgeYbyZFunc(zNewSet(end), 'up'), 'linear', 'extrap');
    [yNewSet(end), zNewSet(end)] = getEdgeYbyZFunc(zNewSet(end), 'up', yNewSet(end));

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Branch 9');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
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

    if ~isreal(yNewSet) || ~isreal(zNewSet)
        warning('检测到复数坐标：yNewSet 或 zNewSet 包含复数');
        plotCrackCoordinates(yNewSet, zNewSet, 'Default Branch');
        % fprintf('--- 原始数据 ---\n');
        % origy
        % origz
        % fprintf('--- 处理后数据 ---\n');
        % yNewSet
        % zNewSet
    end
end

end



