%% Demo script to compare two sets of crack growth parameters (Units: Meters)
clear; clc;

% Baseline parameters from validate_pred_a_N.m
theta1 = [-12.5479, 3.9835, 0.33274645, 6.674963682];

% Modified parameters (hypothetical comparison)
% theta2 = [-12.6673, 4.2264, 0.1142, 6.6750];
theta2 = [-12.2554, 3.2548, 0.9724, 6.6750];
labels = {'Validated Parameters', 'Alternative Estimation'};

% Call the plotting function (Range adjusted for Meters: 5 to 100 MPa*m^0.5)
% Note: This calls the function defined in plot_dadn_comparison.m
plot_dadn_comparison(theta1, theta2, [5, 40], 0.1, labels);

fprintf('Plot generated successfully.\n');
