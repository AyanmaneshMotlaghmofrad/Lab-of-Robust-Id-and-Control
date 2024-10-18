% Create figure
figure('Color', 'white');
hold on;

% Define theta for the parametric curve
theta = linspace(0, 2 * pi, 500);
r = 2 + sin(3 * theta);  % Radius function for a non-convex shape

% Define the non-convex curvy shape
x = r .* cos(theta) + 3;  % Center around (3, 3)
y = r .* sin(theta) + 3;

% Plot the non-convex region using RGB colors
fill(x, y, [0.678 0.847 0.902], 'EdgeColor', [0 0 1], 'LineWidth', 2, 'FaceAlpha', 0.7); % Light blue

% Set axis labels and limits
xlabel('\theta_1', 'FontSize', 14);
ylabel('\theta_2', 'FontSize', 14);
xlim([0 6]);
ylim([0 6]);

% Add grid, title, and aspect ratio
grid on;
title('Non-Convex Graphical Representation of \(\mathbb{D}_\theta\)', 'FontSize', 16);
axis equal;

% Hold off to finish the plot
hold off;
