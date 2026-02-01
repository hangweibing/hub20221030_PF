function plot_dadn_comparison_3d(theta1, theta2, dk_range, r_range, labels)
% PLOT_DADN_COMPARISON_3D Compare two sets of parameters in 3D
% xy-axes: R (Stress Ratio) and Delta K (Log scale)
% z-axis: da/dN (Crack Growth Rate, Log scale)

% Default parameters
if nargin < 3 || isempty(dk_range), dk_range = [5, 40]; end
if nargin < 4 || isempty(r_range), r_range = [-1, 0.9]; end
if nargin < 5 || isempty(labels), labels = {'Set 1', 'Set 2'}; end

% Create Grid
r_vec = linspace(r_range(1), r_range(2), 50);
dk_vec = logspace(log10(dk_range(1)), log10(dk_range(2)), 50);
[R_mesh, DK_mesh] = meshgrid(r_vec, dk_vec);

% Relationship: Kmax = DeltaK / (1 - R)
KMAX_mesh = DK_mesh ./ (1 - R_mesh);

% Model: da/dN = 10^theta(1) * dk^theta(2) * max(kmax/theta(4) - 1, 0)^theta(3)
calc_dadn = @(t, dk, kmax) (10^t(1)) * (dk.^t(2)) .* (max(kmax./t(4) - 1, 0).^t(3));

dadn1 = calc_dadn(theta1, DK_mesh, KMAX_mesh);
dadn2 = calc_dadn(theta2, DK_mesh, KMAX_mesh);

% --- Simple Visualization ---
figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

% Surface 1: Blue
s1 = surf(R_mesh, DK_mesh, dadn1, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', labels{1});

% Surface 2: Red
s2 = surf(R_mesh, DK_mesh, dadn2, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', labels{2});

% Axes Scaling
set(gca, 'YScale', 'log', 'ZScale', 'log');
view(45, 30);
grid on;

% Labels
xlabel('Stress Ratio R');
ylabel('\DeltaK (MPa \cdot m^{0.5})');
zlabel('da/dN (m/cycle)');
title('3D Comparison: da/dN vs R vs \DeltaK');

legend('Location', 'northeast');
hold off;
end
