function da_dN = crack_model_PSR(X, deltaK, R)
% X: N × 4 参数矩阵
% [log(theta1), theta2, theta3, k2]

log_theta1 = X(:,1);
theta2     = X(:,2);
theta3     = X(:,3);
k2         = X(:,4);

% 参数转换
theta1 = 10.^log_theta1;

% Kmax
Kmax = deltaK ./ (1 - R);

% bracket term
bracket_term = Kmax ./ k2 - 1;
bracket_term = max(bracket_term, 0);

% 裂纹扩展速率
da_dN = theta1 .* (deltaK.^theta2) .* (bracket_term.^theta3);

% 防止 NaN / Inf
da_dN(~isfinite(da_dN)) = 0;
end