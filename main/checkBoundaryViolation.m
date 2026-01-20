function [hasViolation, violationInfo] = checkBoundaryViolation(ySet, zSet)
%CHECKBOUNDARYVIOLATION 检查裂纹坐标是否超出物理边界
%
% 输入参数：
%   ySet - y坐标数组
%   zSet - z坐标数组
%
% 输出参数：
%   hasViolation - 布尔值，true表示有坐标点超出边界
%   violationInfo - 违规信息结构体数组，包含违规的坐标点信息
%
% 边界定义：
%   区域1: z∈[30,35]，上边界y=13，下边界y=7
%   区域2: z∈[35,37.82842712]，圆弧边界
%   区域3: z∈[37.82842712,65]，上边界y=11，下边界y=9

    hasViolation = false;
    violationInfo = struct('point', {}, 'z', {}, 'y', {}, 'region', {}, 'reason', {});

    % 设置浮点精度余量，避免由于计算精度导致的误报
    tolerance = 1e-3;

    % 圆心坐标和半径（与plotCrackCoordinates.m保持一致）
    downCenter = [37.82842712, 6];      % 下圆角圆心坐标
    upCenter = [37.82842712, 14];       % 上圆角圆心坐标
    radius = 3;                         % 圆弧半径

    % 检查每个坐标点
    for i = 1:length(ySet)
        y = ySet(i);
        z = zSet(i);

        % 根据z坐标确定所在区域
        if z >= 30 && z <= 35
            % 区域1: 两条平行线
            if y > 13 + tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域1', 'reason', sprintf('超出上边界y=13 (实际y=%.6f)', y));
            elseif y < 7 - tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域1', 'reason', sprintf('超出下边界y=7 (实际y=%.6f)', y));
            end

        elseif z > 35 && z < 37.82842712
            % 区域2: 圆弧区域
            % 计算该z坐标对应的上边界和下边界y值
            y_upper_boundary = upCenter(2) - sqrt(radius^2 - (z - upCenter(1))^2);
            y_lower_boundary = downCenter(2) + sqrt(radius^2 - (z - downCenter(1))^2);

            % 检查点是否在边界内部（上边界和下边界之间），考虑浮点精度余量
            if y > y_upper_boundary + tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域2', 'reason', sprintf('超出上边界y=%.6f (实际y=%.6f)', y_upper_boundary, y));
            elseif y < y_lower_boundary - tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域2', 'reason', sprintf('超出下边界y=%.6f (实际y=%.6f)', y_lower_boundary, y));
            end

        elseif z >= 37.82842712 && z <= 65
            % 区域3: 两条平行线
            if y > 11 + tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域3', 'reason', sprintf('超出上边界y=11 (实际y=%.6f)', y));
            elseif y < 9 - tolerance
                hasViolation = true;
                violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '区域3', 'reason', sprintf('超出下边界y=9 (实际y=%.6f)', y));
            end

        else
            % z坐标超出整个边界范围
            hasViolation = true;
            violationInfo(end+1) = struct('point', i, 'z', z, 'y', y, 'region', '超出边界范围', 'reason', 'z坐标超出[30,65]范围');
        end
    end
end