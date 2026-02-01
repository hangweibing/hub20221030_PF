%% 3D Comparison Demo
clear; clc;

% True params
theta_true = [-12.5479, 3.9835, 0.3327, 6.6749];

% Test params
theta_test = [-12.2554, 3.2548, 0.9724, 6.6750];

labels = {'True Parameters', 'Test Parameters'};

% Plot in 3D: R from -1 to 0.9, DeltaK from 5 to 40
plot_dadn_comparison_3d(theta_true, theta_test, [5, 40], [-1, 0.9], labels);

fprintf('3D Comparison Plot generated.\n');
