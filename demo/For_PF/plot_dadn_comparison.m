function plot_dadn_comparison(theta1, theta2, dk_range, R, labels)
    % PLOT_DADN_COMPARISON Compare two sets of fatigue crack growth parameters
    %
    % Syntax:
    %   plot_dadn_comparison(theta1, theta2)
    %   plot_dadn_comparison(theta1, theta2, dk_range, R, labels)
    %
    % Input:
    %   theta1   - [log10_C, n, p, k] for first set
    %   theta2   - [log10_C, n, p, k] for second set
    %   dk_range - [min, max] range for Delta K (optional)
    %   R        - Load ratio (default 0.1)
    %   labels   - Cell array of strings for legend (optional)

    % Default parameters
    if nargin < 3 || isempty(dk_range)
        dk_range = [5, 100]; % Default range for Delta K (MPa*sqrt(m))
    end
    if nargin < 4 || isempty(R)
        R = 0.1; % Default load ratio
    end
    if nargin < 5 || isempty(labels)
        labels = {'Parameter Set 1', 'Parameter Set 2'};
    end

    % Generate Delta K values (logarithmic spacing)
    dk_vec = logspace(log10(dk_range(1)), log10(dk_range(2)), 200);
    
    % Calculate Kmax based on Delta K and Load Ratio R
    % Delta K = Kmax - Kmin, Kmin = R * Kmax => Delta K = Kmax(1 - R)
    kmax_vec = dk_vec / (1 - R);

    % Model: da/dN = 10^theta(1) * dk^theta(2) * max(kmax/theta(4) - 1, 0)^theta(3)
    calc_dadn = @(t, dk, kmax) (10^t(1)) * (dk.^t(2)) .* (max(kmax./t(4) - 1, 0).^t(3));

    dadn1 = calc_dadn(theta1, dk_vec, kmax_vec);
    dadn2 = calc_dadn(theta2, dk_vec, kmax_vec);

    % --- Premium Visualization ---
    fig = figure('Color', 'w', 'Position', [100, 100, 900, 650]);
    hold on;
    
    % Colors (Premium Dark/Modern Palette)
    color1 = [0.00, 0.45, 0.74]; % Modern Blue
    color2 = [0.85, 0.33, 0.10]; % Modern Orange
    gridColor = [0.94, 0.94, 0.94];

    % Plotting with rich lines
    p1 = semilogx(dk_vec, dadn1, 'LineWidth', 2.5, 'Color', color1, 'DisplayName', labels{1});
    p2 = semilogx(dk_vec, dadn2, 'LineWidth', 2.5, 'Color', color2, 'DisplayName', labels{2});

    % Improve appearance
    set(gca, 'XScale', 'log', 'YScale', 'linear');
    grid on;
    set(gca, 'GridColor', [0.4 0.4 0.4], 'GridAlpha', 0.15, 'MinorGridAlpha', 0.1);
    
    % Labels and Title
    xlabel('\DeltaK (MPa \cdot m^{0.5})', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('da/dN (m/cycle)', 'FontSize', 14, 'FontWeight', 'bold');
    title_str = sprintf('Comparison of Crack Growth Models (R = %.2f)', R);
    title(title_str, 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);

    % Legend
    lgd = legend('Location', 'best');
    set(lgd, 'FontSize', 12, 'Box', 'on', 'EdgeColor', [0.8 0.8 0.8]);

    % Parameter Info Box (Optional but adds value)
    txt1 = sprintf('P1: [%.2f, %.2f, %.2f, %.2f]', theta1(1), theta1(2), theta1(3), theta1(4));
    txt2 = sprintf('P2: [%.2f, %.2f, %.2f, %.2f]', theta2(1), theta2(2), theta2(3), theta2(4));
    annotation('textbox', [0.15, 0.15, 0.3, 0.1], 'String', {['\color[rgb]{0.00,0.45,0.74}', txt1], ['\color[rgb]{0.85,0.33,0.10}', txt2]}, ...
        'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5, 'Margin', 5);

    % Better ticks
    set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out', 'TickLength', [0.02, 0.02]);
    
    hold off;
end
