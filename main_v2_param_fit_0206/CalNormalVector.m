%% ===================================================================
%% 函数名称：CalNormalVector
%% 功能描述：基于裂纹几何形状计算法向量
%% ===================================================================
function [normalNormalizeVectorSet] = CalNormalVector(yRegSet, zRegSet, ksiRegSet)

% =========================================================================
% 输入参数：
%   yRegSet:         裂纹轮廓的y坐标集 (21个节点)
%   zRegSet:         裂纹轮廓的z坐标集 (21个节点)
%   ksiRegSet:       参数化坐标系的参数值 (21个节点)
%
% 输出参数：
%   normalNormalizeVectorSet: 归一化的法向量集合 (2×21矩阵)
% =========================================================================

%% 三次样条插值创建平滑曲线
% 使用二阶连续的三次样条插值创建y坐标的平滑曲线
csy = csape(ksiRegSet, yRegSet, 'second');

% 使用二阶连续的三次样条插值创建z坐标的平滑曲线
csz = csape(ksiRegSet, zRegSet, 'second');

%% 数值微分计算切向量
% 设置微分步长，用于数值微分计算
deltan = 0.00001;

% 计算y方向的切向量分量（中心差分公式）
ytangent = (ppval(csy, ksiRegSet+deltan) - ppval(csy, ksiRegSet-deltan)) / (2*deltan);

% 计算z方向的切向量分量（中心差分公式）
ztangent = (ppval(csz, ksiRegSet+deltan) - ppval(csz, ksiRegSet-deltan)) / (2*deltan);

% 注意：另一种计算方式（已注释）
% ytangent = (ppval(csy, ksiRegSet+0.0001) - ppval(csy, ksiRegSet-0.0001)) / 0.0002;

%% 计算法向量
% 组合切向量（添加负号调整方向）
tangentVector = -[ytangent; ztangent];

% 归一化切向量（单位化）
tangentNormalizeVectorSet = tangentVector ./ (sqrt(tangentVector(1,:).^2 + tangentVector(2,:).^2));

% 计算法向量（垂直于切向量，逆时针旋转90度）
normalNormalizeVectorSet = [tangentNormalizeVectorSet(2,:); -tangentNormalizeVectorSet(1,:)];

end
