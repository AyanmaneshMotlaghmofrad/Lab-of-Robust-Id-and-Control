% Parameters for D_theta (an ellipse)
center = [0.5, 0.5];  % Center of the ellipse
width = 1.0;          % Width of the ellipse
height = 0.5;         % Height of the ellipse
rotation_angle = 30;  % Rotation angle in degrees

% Generate points for the ellipse
theta = linspace(0, 2 * pi, 100);
x = (width / 2) * cos(theta);
y = (height / 2) * sin(theta);

% Rotation matrix
rotation_matrix = [cosd(rotation_angle), -sind(rotation_angle);
                   sind(rotation_angle),  cosd(rotation_angle)];

% Rotate the ellipse
rotated_points = rotation_matrix * [x; y];
rotated_x = rotated_points(1, :) + center(1);
rotated_y = rotated_points(2, :) + center(2);

% Calculate the minimum and maximum values
min_theta1 = min(rotated_x);
max_theta1 = max(rotated_x);
min_theta2 = min(rotated_y);
max_theta2 = max(rotated_y);

% Create a figure
figure;

% Plot the rotated ellipse representing D_theta
plot(rotated_x, rotated_y, 'b', 'DisplayName', '$\mathbb{D}_\theta$');
hold on;

% Create a rectangle representing the bounding box
bounding_box = rectangle('Position', [min_theta1, min_theta2, max_theta1 - min_theta1, max_theta2 - min_theta2], ...
                         'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'none');

% Set labels for axes with minimum and maximum values
xlabel('$\theta_1$', 'Interpreter', 'latex', 'FontSize', 14);  % Adjust font size
ylabel('$\theta_2$', 'Interpreter', 'latex', 'FontSize', 14);  % Adjust font size
title('Bounding box and $\mathbb{D}_\theta$', 'Interpreter', 'latex', 'FontSize', 16);  % Adjust font size

% Annotate minimum and maximum points with underlining and overlining
function annotate_with_line(x, y, text_str, line_y_offset, color)
    % Annotate the text
    hText = text(x, y, text_str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Get the text extent
    ext = get(hText, 'Extent');
    
    % Calculate the line width based on the text width
    line_width = ext(3);  % Width of the text
    
    % Draw a line for underlining or overlining
    line([x - line_width / 2, x + line_width / 2], ...
         [y + line_y_offset, y + line_y_offset], ...
         'Color', color, 'LineWidth', 1);
end

% Underline and overline annotations
annotate_with_line(min_theta1, center(2), '$\theta_{1}$', -0.05, 'k');
annotate_with_line(max_theta1, center(2), '$\theta_{1}$', -0.05, 'k');
annotate_with_line(center(1), min_theta2, '$\theta_{2}$', 0.05, 'k');
annotate_with_line(center(1), max_theta2, '$\theta_{2}$', 0.05, 'k');

% Show grid and set aspect ratio to equal
grid on;
axis equal;
xlim([-0.5, 1.5]);
ylim([-0.5, 1.5]);
hold off;
