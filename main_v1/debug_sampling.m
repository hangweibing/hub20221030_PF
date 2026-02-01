function [log_theta1_, theta2, theta3, k2] = debug_sampling()
%% DEBUG_SAMPLING 用于调试模式的参数采样函数
%
% 该函数直接从 'AM-TC4-GRO.xlsx' 文件中随机抽取一组现有的参数组合
% 用于调试模式，避免 Copula 采样可能带来的不确定性
%
% 输出参数:
%   log_theta1_ - logD 参数 (对应训练数据第一列)
%   theta2      - A 参数 (对应训练数据第二列)
%   theta3      - delta_kthr 参数 (对应训练数据第三列)
%   k2          - p 参数 (对应训练数据第四列)

    % 读取训练数据
    filename = 'AM-TC4-GRO.xlsx';
    try
        T = readtable(filename);
        X_train = table2array(T);
        [N_train, d] = size(X_train);

        if d < 4
            error('训练数据维度不足，需要至少4列参数');
        end

    catch ME
        error('无法读取训练数据文件 %s: %s', filename, ME.message);
    end

    % 从现有数据中随机抽取一行
    random_idx = randi(N_train);

    % 提取参数
    log_theta1_ = X_train(random_idx, 1);
    theta2 = X_train(random_idx, 2);
    theta3 = X_train(random_idx, 3);
    k2 = X_train(random_idx, 4);

    
end