function K = sim_K_func(a, P, W, B)
% SIM_K_FUNC C(T) Specimen Stress Intensity Factor calculation
% 依据 ASTM E647 A1.5.1.1 公式实现
% a: 裂纹长度 crack length (mm)
% P: 载荷 load (N) 或 载荷幅值 delta P
% W: 试样宽度 width (mm)
% B: 试样厚度 thickness (mm)

% 1. 计算无量纲裂纹长度 alpha
alpha = a / W;

% % 校验有效范围: 根据图中说明 alpha >= 0.2
% if alpha < 0.2
%     warning('裂纹长度比例 a/W < 0.2，超出该经验公式的严谨适用范围。');
% elseif alpha >= 1
%     K = NaN;  % 裂纹长度超过试样宽度，物理上不可能，返回 NaN
%     return;
% end

% 2. 计算多项式部分 (Polynomial expression)
f_alpha = (0.886 + 4.64*alpha - 13.32*alpha^2 + 14.72*alpha^3 - 5.6*alpha^4);

% 3. 计算前置几何因子部分
geometry_term = (2 + alpha) / (1 - alpha)^(1.5);

% 4. 完整公式计算: K = (P / (B * sqrt(W))) * geometry_term * f_alpha
% 注意：单位输出为 MPa*sqrt(mm)
K = (P / (B * sqrt(W))) * geometry_term * f_alpha;
end