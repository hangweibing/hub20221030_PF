%% ===================================================================
%% 函数名称：crackRegular5Func
%% 功能描述：生成平滑规则的裂纹轮廓（等弧长分布）
%% ===================================================================
function [yRegSet, zRegSet, ksiIniSet, splitted] = crackRegular5Func(yIniSet, zIniSet, nRegPoint, indexNeedIniNatCoord)

% =========================================================================
% 输入参数：
%   yIniSet:             初始裂纹y坐标集
%   zIniSet:             初始裂纹z坐标集
%   nRegPoint:           目标节点数量（通常为21）
%   indexNeedIniNatCoord: 是否需要初始自然坐标
%
% 输出参数：
%   yRegSet:             等弧长分布的y坐标集
%   zRegSet:             等弧长分布的z坐标集
%   ksiIniSet:           初始自然坐标参数（如果需要）
%   splitted:             裂纹是否发生分裂
%
% 功能说明：
%   将任意形状的裂纹轮廓转换为平滑的等弧长分布曲线
%   确保几何连续性和数值稳定性
% =========================================================================

%% 初始化参数
% 工件几何约束参数
downCenter = [6, 37.82842712];     % 下圆角圆心
r = 3;                             % 圆角半径
downCenter1 = [6, 37.82842712];    % 下圆角圆心（备用）
%% ===================================================================
%% 第一阶段：前缘分裂处理
%% ===================================================================

% 处理裂纹前缘与圆角的相交情况
[yIniSet, zIniSet, splitted] = resplit_front(yIniSet, zIniSet, downCenter1, r);

%% ===================================================================
%% 第二阶段：智能节点选择
%% ===================================================================

% 获取原始节点数量
nOriPoint = size(yIniSet, 2);

% 为智能采样选择参考节点位置
% 避免选择与圆角相交的区域，确保几何稳定性
mid_left_point_index = ceil(nOriPoint/3);    % 左侧中间点索引
mid_right_point_index = floor(2*nOriPoint/3); % 右侧中间点索引

%% 检查节点与圆角的相交情况
judge = zeros(1, nOriPoint);

% 判断每个节点是否在圆角区域内
for i = 1:length(yIniSet)
    judge(i) = isPointInCircleFunc([yIniSet(i), zIniSet(i)], downCenter, r);
end

% 获取与圆角相交的节点索引
intersected_coords = find(judge == true);
%% 根据相交情况选择节点
if isempty(intersected_coords)
    %% 子情况1：无相交节点，标准四点采样
    % 选择起始点、左侧中间点、右侧中间点、结束点
    try
        zIniSet = zIniSet([1, mid_left_point_index, mid_right_point_index, end]);
        yIniSet = yIniSet([1, mid_left_point_index, mid_right_point_index, end]);
    catch
        error()
    end
else
    %% 子情况2：存在相交节点，避开圆角区域
    % 获取相交区域的边界
    left_insec_coord = intersected_coords(1);   % 最左侧相交点
    right_insec_coord = intersected_coords(end); % 最右侧相交点

    % 自适应调整中间点位置，避免选择圆角区域
    % 调整中间点位置，避免与圆角区域重叠
    if mid_left_point_index >= left_insec_coord
        mid_left_point_index = left_insec_coord - 1;  % 向左移动到安全区域
    end
    if mid_right_point_index <= right_insec_coord
        mid_right_point_index = right_insec_coord + 1; % 向右移动到安全区域
    end

    %% 根据调整后的中间点进行最终节点选择
    if mid_left_point_index > 1  % 左侧中间点不在边界
        if mid_right_point_index < nOriPoint  % 右侧中间点不在边界
            %% 四点采样：标准情况
            zIniSet = zIniSet([1, mid_left_point_index, mid_right_point_index, end]);
            yIniSet = yIniSet([1, mid_left_point_index, mid_right_point_index, end]);
        else
            %% 三点采样：右侧点在边界附近
            if left_insec_coord <= 4
                % 使用固定前三个点
                zIniSet = zIniSet([1, 2, 3, end]);
                yIniSet = yIniSet([1, 2, 3, end]);
            else
                % 计算两个中间点
                mid_left_point_1_index = floor((2/3 + left_insec_coord/3));
                mid_left_point_2_index = ceil(1/3 + 2*left_insec_coord/3);
                zIniSet = zIniSet([1, mid_left_point_1_index, mid_left_point_2_index, end]);
                yIniSet = yIniSet([1, mid_left_point_1_index, mid_left_point_2_index, end]);
            end
        end
    else  % 左侧中间点在边界
        if right_insec_coord >= nOriPoint-3
            % 使用固定后三个点
            zIniSet = zIniSet([1, nOriPoint-2, nOriPoint-1, end]);
            yIniSet = yIniSet([1, nOriPoint-2, nOriPoint-1, end]);
        else
            % 计算两个中间点（偏向右侧）
            mid_right_point_1_index = floor(2/3*right_insec_coord + nOriPoint/3);
            mid_right_point_2_index = ceil(1/3*right_insec_coord + 2*nOriPoint/3);
            zIniSet = zIniSet([1, mid_right_point_1_index, mid_right_point_2_index, end]);
            yIniSet = yIniSet([1, mid_right_point_1_index, mid_right_point_2_index, end]);
        end
    end
end


%% ===================================================================
%% 第三阶段：B样条曲线重建
%% ===================================================================

% 使用B样条创建平滑曲线
cs = cscvn([zIniSet; yIniSet]);  % B样条曲线拟合

% 调试绘图代码（已注释）
% fnplt(cs)
% legend('原始曲线')
% hold on;

%% ===================================================================
%% 第四阶段：等弧长参数化
%% ===================================================================

% 已注释的等参数方法（备用）
% equi_para = linspace(0, cs.breaks(end), nRegPoint);
% test = fnval(cs, equi_para);
% plot(test(1,:), test(2,:), '*')
% legend('原始曲线')

%% 创建密集的曲线采样点，用于计算弧长
n_dense_point = 3000;  % 密集采样点数量
dense_para = linspace(0, cs.breaks(end), n_dense_point);
dense_points = fnval(cs, dense_para);  % 计算密集采样点的坐标

% 调试绘图代码（已注释）
% plot(dense_points(1,:), dense_points(2,:), '^')

%% 计算总弧长
s_tot = 0;
for i = 1:n_dense_point-1
    % 计算相邻两点间的距离并累加
    s_tot = s_tot + sqrt((dense_points(:,i+1) - dense_points(:,i))' * (dense_points(:,i+1) - dense_points(:,i)));
end

%% 创建等弧长分布的间隔点
s_interval = linspace(0, s_tot, nRegPoint);  % 等弧长间隔

%% 初始化输出数组
para_ori = cs.breaks;  % B样条参数节点
yRegSet = zeros(1, nRegPoint);  % 等弧长y坐标
zRegSet = zeros(1, nRegPoint);  % 等弧长z坐标

%% 设置边界点（起始点和结束点）
zRegSet(1) = dense_points(1, 1);    % 起始点z坐标
yRegSet(1) = dense_points(2, 1);    % 起始点y坐标
zRegSet(end) = dense_points(1, end); % 结束点z坐标
yRegSet(end) = dense_points(2, end); % 结束点y坐标

%% 初始化弧长计算变量
s_ori = zeros(1, nOriPoint);  % 原始节点的弧长位置
s_ori(1) = 0;                 % 起始点弧长为0
s_ori(end) = s_tot;           % 结束点弧长为总弧长

%% 循环变量初始化
s_cur = 0;      % 当前累积弧长
s_his = 0;      % 历史累积弧长
cur_s_index = 2;  % 当前搜索等弧长点的索引
cur_p_index = 2;  % 当前搜索原始参数点的索引
%% 主循环：计算等弧长分布的节点
for i = 1:n_dense_point-1
    % 更新弧长累积
    s_his = s_cur;  % 保存上一段的弧长
    s_cur = s_cur + sqrt((dense_points(:,i+1) - dense_points(:,i))' * (dense_points(:,i+1) - dense_points(:,i))); % 累加当前段弧长

    %% 检查是否到达等弧长点
    if s_cur > s_interval(cur_s_index) && s_his < s_interval(cur_s_index)
        % 在当前段内插值找到等弧长点
        zRegSet(cur_s_index) = (dense_points(1,i) + dense_points(1,i+1)) / 2;  % z坐标插值
        yRegSet(cur_s_index) = (dense_points(2,i) + dense_points(2,i+1)) / 2;  % y坐标插值
        cur_s_index = cur_s_index + 1;  % 移动到下一个等弧长点
    end

    %% 可选：计算原始节点的自然坐标
    if indexNeedIniNatCoord
        % 将原始B样条参数转换为弧长参数
        if para_ori(cur_p_index) >= dense_para(i) && para_ori(cur_p_index) < dense_para(i+1)
            % 计算原始节点在弧长参数下的位置
            s_ori(cur_p_index) = s_his + ((fnval(cs, para_ori(cur_p_index)) - dense_points(:,i))' * ...
                                 (fnval(cs, para_ori(cur_p_index)) - dense_points(:,i)));
            cur_p_index = cur_p_index + 1;
        end
    end
end

%% 调试绘图（已注释）
% plot(zRegSet, yRegSet, '^')

%% 输出自然坐标参数（如果需要）
if indexNeedIniNatCoord
    ksiIniSet = (s_ori ./ s_tot);  % 归一化到[0,1]区间
end

%% 函数结束
% 注意：y和z坐标已按等弧长分布重新排列，确保几何连续性

end

