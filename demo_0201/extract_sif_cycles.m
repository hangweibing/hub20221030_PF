%% 应力循环提取与应力强度因子计算程序
clear; clc;

% (1) & (2) 参数设置
num_cycles = 100; % 手动设置需要提取几个循环
ref_sigma = 100;  % 基准应力 (MPa)
ref_K = 15;       % 基准应力对应的 K

% 加载数据 (确保 AsteixSpectraData.mat 在当前路径或搜索路径中)c
load('AsteixSpectraData.mat', 'spectra'); 

% 预处理载荷谱 (参考原程序逻辑：重复并去除首点)
spectra = repmat(reshape(spectra, 1, []), 1, 7);
spectra = spectra(2:end);

% (3) & (4) 计算应力强度因子并组成向量
sif_factor = ref_K / ref_sigma;
sif_results = zeros(num_cycles, 2);

for i = 1:num_cycles
    Smin = spectra(2*i-1);
    Smax = spectra(2*i);
    
    Kmax = Smax * sif_factor;
    deltaK = (Smax - Smin) * sif_factor;
    
    sif_results(i, :) = [deltaK, Kmax];
end

% (5) 保存结果
save('sif_cycles_extracted.mat', 'sif_results');

fprintf('提取完成，共 %d 个循环，结果已保存至 sif_cycles_extracted.mat\n', num_cycles);
disp('前 5 行结果 [deltaK, Kmax]:');
disp(sif_results(1:min(5, num_cycles), :));
