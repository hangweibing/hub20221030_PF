

clc
clear all
close all

load('AsteixSpectraData.mat')

% 输入：spectra（长度已保证为偶数）
spectra = spectra(:);

thr = 30;                 % 你的阈值
n = floor(numel(spectra)/2)*2;
s = spectra(1:n);

m = (s(1:2:end) + s(2:2:end))/2;   % 每个循环均值
dm = abs(diff(m));                % 与上一循环均值差

idxCycle = find(dm > thr) + 1;     % 第 idxCycle 个循环发生突变（相对前一个）
A = 2*idxCycle - 1;                % 映射回 spectra 的位置：该循环第一个点的下标




spectra = spectra(:);

% ---- 参数：量化步长（按你的单位设置）----
mean_step = 1;     % 均值按 1（MPa等）分辨率统计；需要更细就改 0.1
R_step    = 0.01;  % R 按 0.01 分辨率统计

% ---- 两点一循环 ----
n = floor(numel(spectra)/2)*2;
s = spectra(1:n);
hi = s(1:2:end);
lo = s(2:2:end);

m = (hi + lo)/2;                 % 循环均值
R = lo ./ hi;                    % 应力比：σmin/σmax（这里假设第1点是max，第2点是min）

% 若你不确定第1点一定是max，可用下面这两行更稳妥：
% hi = max([s(1:2:end), s(2:2:end)], [], 2);
% lo = min([s(1:2:end), s(2:2:end)], [], 2);
% m  = (hi + lo)/2;  R = lo./hi;

% ---- 处理异常：hi=0 导致 R=Inf/NaN ----
bad = ~isfinite(R);
R(bad) = NaN;

% ---- 量化（避免浮点噪声导致“种类爆炸”）----
mq = round(m/mean_step)*mean_step;
Rq = round(R/R_step)*R_step;

% ---- 统计均值种类与频率 ----
[uM,~,iM] = unique(mq);
cntM = accumarray(iM, 1);
freqM = cntM / numel(mq);
T_mean = table(uM, cntM, freqM, 'VariableNames', {'Mean', 'Count', 'Freq'});
T_mean = sortrows(T_mean, 'Count', 'descend');

% ---- 统计应力比种类与频率（忽略 NaN）----
maskR = ~isnan(Rq);
[uR,~,iR] = unique(Rq(maskR));
cntR = accumarray(iR, 1);
freqR = cntR / sum(maskR);
T_R = table(uR, cntR, freqR, 'VariableNames', {'R', 'Count', 'Freq'});
T_R = sortrows(T_R, 'Count', 'descend');

% ---- 输出 ----
nMeanTypes = height(T_mean);
nRTypes    = height(T_R);

disp(nMeanTypes);  % 均值种类数
disp(T_mean);      % 均值-次数-频率表

disp(nRTypes);     % 应力比种类数
disp(T_R);         % R-次数-频率表


spectra = spectra(:);

% ===== 参数（你只需要调这两个）=====
low_thr = 65;     % 低载荷阈值（MPa，按你的图 60~70 都合理）
min_len = 5;      % 至少连续多少个 cycle 才算一个“段”

% ===== 1. 标记低载荷点 =====
isLow = spectra < low_thr;

% ===== 2. 找连续区间 =====
d = diff([0; isLow; 0]);
seg_start = find(d == 1);
seg_end   = find(d == -1) - 1;

seg_len = seg_end - seg_start + 1;

% ===== 3. 过滤太短的段 =====
valid = seg_len >= min_len;

seg_start = seg_start(valid);
seg_end   = seg_end(valid);
seg_len   = seg_len(valid);

% ===== 4. 结果 =====
nSegment = numel(seg_start);   % 这样的段的个数

disp(['识别到的低载荷段数量 = ', num2str(nSegment)]);



