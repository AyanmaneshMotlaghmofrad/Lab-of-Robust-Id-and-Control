% Sample data for the input signal u(k)
k = 0:10; % Time indices from 0 to 10
u_k = [1, 2, 3, 5, 4, 6, 5, 7, 6, 8, 7]; % Input signal values

% Sample data for the output signal y(k)
y_k = sin(u_k) + 0.1 * randn(size(u_k)); % Non-linear transformation with noise
delta_eta = 1; % Noise level for output
y_k_upper = y_k + delta_eta; % Upper bound
y_k_lower = y_k - delta_eta; % Lower bound

% Create the figure and axes
figure;

% Input Signal Plot
subplot(2, 1, 1);
plot(k, u_k, 'o-', 'Color', [0.39, 0.58, 0.93], 'LineWidth', 2, 'MarkerSize', 6); % Cornflower Blue
title('Input Signal u(k)', 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('k (Time)', 'FontSize', 12, 'Interpreter', 'latex');
ylabel('u(k)', 'FontSize', 12, 'Interpreter', 'latex');
grid on;

% Output Signal Plot
subplot(2, 1, 2);
hold on;
plot(k, y_k, 'o-', 'Color', [1, 0.39, 0.28], 'LineWidth', 2, 'MarkerSize', 6); % Tomato
plot(k, y_k_upper, '--', 'Color', [0.24, 0.70, 0.44], 'LineWidth', 1.5); % Medium Sea Green
plot(k, y_k_lower, '--', 'Color', [0.24, 0.70, 0.44], 'LineWidth', 1.5); % Medium Sea Green
fill([k fliplr(k)], [y_k_lower fliplr(y_k_upper)], [0.24, 0.70, 0.44], 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Fill between bounds
title('Output Signal y(k) with Bounds', 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$\tilde{y}(k)$', 'FontSize', 12, 'Interpreter', 'latex');  % Changed the horizontal label
ylabel('$\tilde{y}(k)$', 'FontSize', 12, 'Interpreter', 'latex');
grid on;
legend('$\tilde{y}(k)$', '$\tilde{y}(k) + \Delta \eta$', '$\tilde{y}(k) - \Delta \eta$', 'Bounds', 'Location', 'Best', 'Interpreter', 'latex');

% Adjust layout
tight_layout();
