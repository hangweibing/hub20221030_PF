clc; clear; close all;

%% ================== 1. 读取参数文件 ==================
filename = 'NASGRO(H-S).xlsx';  
T = readtable(filename);
param_names = T.Properties.VariableNames;
Xdata = table2array(T);   % N × 4
d = size(Xdata, 2);       % 参数维数

%% ================== 2. 定义参数范围 ==================
lb = min(Xdata, [], 1);
ub = max(Xdata, [], 1);

%% ================== 3. Sobol 样本收敛性分析设置 ==================
sample_list = round(linspace(500, 50000, 10)); % 10 个样本数量
S_all  = zeros(d, numel(sample_list));
ST_all = zeros(d, numel(sample_list));

deltaK = 6;  
R_ratio = 0.5;
Kmax = deltaK / (1 - R_ratio);

%% ================== 4. 循环计算不同样本数的 Sobol 指数 ==================
for k = 1:numel(sample_list)
    N = sample_list(k);
    
    % Sobol 序列采样
    pset = sobolset(2*d, 'Skip', 1000, 'Leap', 100);
    pset = scramble(pset, 'MatousekAffineOwen');
    R_sobol = net(pset, N); % [N x 2d]
    
    VarMin = [lb, lb]; 
    VarMax = [ub, ub];
    R_scaled = VarMin + R_sobol .* (VarMax - VarMin);
    
    A = R_scaled(:, 1:d);
    B = R_scaled(:, d+1:end);
    
    % 模型输出
    YA = NASGRO_model(A, deltaK, Kmax);
    YB = NASGRO_model(B, deltaK, Kmax);
    
    VY = var([YA; YB], 1);
    S  = zeros(d,1);
    ST = zeros(d,1);
    
    for i = 1:d
        AB_i = A;
        AB_i(:, i) = B(:, i);
        YAB_i = NASGRO_model(AB_i, deltaK, Kmax);
        
        S(i)  = mean(YB .* (YAB_i - YA)) / VY;
        ST(i) = 0.5 * mean((YA - YAB_i).^2) / VY;
    end
    
    S_all(:,k)  = S;
    ST_all(:,k) = ST;
end

%% ================== 5. 输出最终结果 (最大样本数) ==================
disp('--- Sobol 敏感性分析结果 (最大样本数) ---');
SobolTable = table(param_names', S_all(:,end), ST_all(:,end), ...
    'VariableNames', {'Parameter','FirstOrder','TotalOrder'});
disp(SobolTable);

%% ================== 6. 可视化原有 Sobol 条形图 ==================
figure('Color','w');
subplot(1,2,1);
bar(S_all(:,end), 'FaceColor',[0.3 0.6 0.9]);
set(gca,'XTick',1:d,'XTickLabel',param_names);
ylabel('First-order Index (S_1)');
title('First-order Sensitivity');
grid on;

subplot(1,2,2);
bar(ST_all(:,end), 'FaceColor',[0.9 0.4 0.4]);
set(gca,'XTick',1:d,'XTickLabel',param_names);
ylabel('Total-order Index (S_T)');
title('Total-order Sensitivity');
grid on;
sgtitle(['Sobol Sensitivity Analysis (N = ', num2str(sample_list(end)), ')']);

%% ================== 7. 样本数收敛性分析图 ==================
figure('Color','w');
colors = lines(d);            % 四种颜色
markers = {'o','s','^','d'}; % 四种散点符号
legend_handles = gobjects(d,1); % 用于 legend

% 子图 1: S1
subplot(1,2,1); hold on;
for i = 1:d
    plot(sample_list, S_all(i,:), '-','Color',colors(i,:),'LineWidth',1.5); % 仅折线
    h = scatter(sample_list, S_all(i,:), 40, colors(i,:), markers{i}, 'filled'); % 散点
    legend_handles(i) = h; % legend 只使用散点对象
end
xlabel('Sample Number');
ylabel('First-order Index (S_1)');
title('S_1 Convergence with Sample Number');
legend(legend_handles, param_names,'Location','best');
grid on;

% 子图 2: ST
subplot(1,2,2); hold on;
for i = 1:d
    plot(sample_list, ST_all(i,:), '-','Color',colors(i,:),'LineWidth',1.5);
    h = scatter(sample_list, ST_all(i,:), 40, colors(i,:), markers{i}, 'filled');
    legend_handles(i) = h;
end
xlabel('Sample Number');
ylabel('Total-order Index (S_T)');
title('S_T Convergence with Sample Number');
legend(legend_handles, param_names,'Location','best');
grid on;
sgtitle('Sobol Index Convergence Analysis');
    
%% ================== NASGRO 模型函数 (log10 输出) ==================
function da_dN_log = NASGRO_model(X, deltaK, Kmax)
    log10_D = X(:,1); p = X(:,2); K_thr = X(:,3); A = X(:,4);
    D = 10.^log10_D;
    effective_dK = max(deltaK - K_thr, 1e-6);
    denominator = max(1 - Kmax ./ A, 1e-6);
    da_dN = D .* (effective_dK).^p ./ (denominator).^(p/2);
    da_dN_log = log10(da_dN + 1e-12);  
    da_dN_log(~isfinite(da_dN_log)) = 0;
end
