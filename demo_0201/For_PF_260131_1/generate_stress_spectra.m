% generate_stress_spectra.m
% 该程序用于生成多段不同应力比的恒幅载荷谱，替代原有的应力谱文件
% 目标文件: AsteixSpectraData.mat (变量名: spectra, 1*743230)

clear; clc;

% (1) 设置目标长度
target_len = 743230;
% 计算可用周期数 (去除第一个点后)
% 为了保证 [smin, smax] 成对出现，去除首点后的长度应为偶数
total_cycles = floor((target_len - 1) / 2);
actual_len = 1 + 2 * total_cycles; % 实际生成的总长度

% (2) 设定段数和参数
K = 5;

% 定义每段的应力比和平均应力数组 (用户定义的参数)
R_array = [0.7, 0.3, 0.5, 0.4, 0.75];
S_avg_array = [100, 60, 80, 40, 110]; % 修正了之前的 typo

% 分配每段的周期数
cycles_per_seg = floor(total_cycles / K);
seg_cycles = repmat(cycles_per_seg, 1, K);
seg_cycles(end) = total_cycles - sum(seg_cycles(1:end-1));

% +1 for the transition point at the end
spectra = zeros(1, actual_len + 1);


% 索引 1: 预设为第一段的 Smax (占位符，在 pred_a_N 中通过 spectra(2:end) 丢弃)
S_max1 = (2 * S_avg_array(1)) / (1 + R_array(1));
spectra(1) = S_max1;

% (3) 循环生成每段应力谱
current_idx = 2;
fprintf('=== 恒幅应力谱生成程序 (优化版) ===\n');
fprintf('目标总长度: %d, 实际生成长度: %d\n', target_len, actual_len);

for i = 1:K
    R = R_array(i);
    S_avg = S_avg_array(i);

    % 计算峰值和谷值
    S_max = (2 * S_avg) / (1 + R);
    S_min = R * S_max;

    % 每段包含整数个周期，每个周期由 [Smin, Smax] 组成
    % 对应索引模式: [偶数, 奇数, 偶数, 奇数...] -> [Smin, Smax, Smin, Smax...]
    num_pts = 2 * seg_cycles(i);

    for j = 1:num_pts
        if mod(current_idx, 2) == 0
            spectra(current_idx) = S_min; % 偶数索引 -> Smin
        else
            spectra(current_idx) = S_max; % 奇数索引 -> Smax
        end
        current_idx = current_idx + 1;
    end

    fprintf('第 %d/%d 段: R=%.2f, S_avg=%.2f, 周期数=%d, 结束索引=%d\n', ...
        i, K, R, S_avg, seg_cycles(i), current_idx-1);
end

% NEW: Append the Smin of the first segment to the end to ensure smooth splicing
% Connector pair will be (spectra(end), spectra(1) of next copy) -> (Smin1, Smax1)
spectra(current_idx) = spectra(2);
fprintf('  -> 添加拼接过渡点: spectra(%d) = %.4f (Loop Smin)\n', current_idx, spectra(current_idx));
actual_len = current_idx; % Update actual length


% (4) 保存文件
save_name = 'AsteixSpectraData_fake.mat';
save(save_name, 'spectra');

% (5) 绘制应力谱折线图
figure;
plot(1:length(spectra), spectra, '-b');
hold on;
% 标记段间界限 (仅示意前100个点或整体)
title('恒幅应力块谱');
grid on;


% ... (existing plotting code) ...
grid on;

% (6) 完整性检查 (模仿 gen_synthetic_data_spectrum.m 的处理方式进行拼接检查)
fprintf('\n=== 载荷谱拼接完整性检查 ===\n');
% 模拟程序中的读取和拼接方式 (取2个周期进行拼接检查连接处)
spectra_long = repmat(spectra, 1, 2);
spectra_check = spectra_long(2:end);

% 提取所有周期的 Smin 和 Smax
% 根据 gen_synthetic_data_spectrum.m:
% Smin_vec = loads_segment(1:2:end);
% Smax_vec = loads_segment(2:2:end);
% 注意：如果总长度为奇数(去除第一个点后)，最后可能会多出一个不成对的点，需处理
len_check = floor(length(spectra_check)/2) * 2;
spectra_check = spectra_check(1:len_check);

Smin_full = spectra_check(1:2:end);
Smax_full = spectra_check(2:2:end);

% 检查 Smax < Smin 的情况
invalid_idx = find(Smax_full < Smin_full);
if ~isempty(invalid_idx)
    warning('Fail: 检测到 %d 个周期的 Smax < Smin！', length(invalid_idx));
    fprintf('  前 5 个异常位置 (索引及值):\n');
    for k = 1:min(5, length(invalid_idx))
        idx = invalid_idx(k);
        fprintf('  Cycle %d: Smin=%.4f, Smax=%.4f (Diff=%.4f)\n', ...
            idx, Smin_full(idx), Smax_full(idx), Smax_full(idx)-Smin_full(idx));
    end
else
    fprintf('Pass: 拼接检查通过，所有周期均满足 Smax >= Smin。\n');
end

% 统计 Smin 和 Smax 的范围

fprintf('生成成功！\n');
fprintf('文件已保存为: %s\n', save_name);


