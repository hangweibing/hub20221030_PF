%% Simplified Merge Model Parameters
% Reads pop_now and updates f1, f2 based on log_theta1 and theta2 formulas.
% Mapping for pop_now{1,1} (Cols 2-6):
% Col 2: k1, Col 3: k2, Col 4: f1, Col 5: f3, Col 6: f2

clear; clc;
data_dir = 'f:\研究生\科研\粒子滤波\hub20221030_PF\demo\For_PF_260125_1';
file_path = fullfile(data_dir, 'AM-TC4-GRO_260123.mat');

if exist(file_path, 'file')
    load(file_path); % Loads pop_now (1x6 cell)
else
    error('File not found: %s', file_path);
end

k1 = pop_now{1,1}(:, 2);
f1 = pop_now{1,1}(:, 4);
f3 = pop_now{1,1}(:, 5);
f2 = pop_now{1,1}(:, 6);

% Update values in-place: 
% Col 4 (f1) -> log_theta1 = f1 - f3 * log10(k1)
% Col 6 (f2) -> theta2 = f2 + f3
pop_now{1,1}(:, 4) = f1 - f3 .* log10(k1); % log_theta1
pop_now{1,1}(:, 6) = f2 + f3;              % theta2
pop_now{1,1}(:, 5) = f3;                   % theta3

p_overall = pop_now{1,6};
if length(p_overall) == 5
    % Index mapping: 1:k1, 2:k2, 3:f1, 4:f3, 5:f2
    f1 = p_overall(3);
    f3 = p_overall(4);
    f2 = p_overall(5);
    p_overall(3) = f1 - f3 * log10(p_overall(1)); % Update f1 to log_theta1
    p_overall(4) = f3;
    p_overall(5) = f2 + f3;                       % Update f2 to theta2
    pop_now{1,6} = p_overall;
end

% Save results back in the same structure
output_path = fullfile(data_dir, 'AM-TC4-GRO_260123_merged.mat');
save(output_path, 'pop_now');
fprintf('Updated pop_now saved to: %s\n', output_path);
