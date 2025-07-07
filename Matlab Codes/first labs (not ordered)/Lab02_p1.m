clear all
close all
clc

%% Problem Definition
z = tf('z',1);
u_tilde = [4, -3, 2, 1]';
theta1 = -0.5;
theta2 = 2;
theta = [theta1, theta2];
e = [0.05, -0.25, 0.3, -0.5]';
delta_e = 0.5;

%% Simulation
t_min = 2;
N = 4;
y_tilde = zeros(4,1);
for t = t_min:N
    y_tilde(t) = -theta1 * y_tilde(t-1) + theta2 * u_tilde(t) + e(t);
end

%% Polytope Plot
theta1_vals = linspace(-3, 3, 1000);
figure(1), hold on;
for t = t_min:N
    theta2_max = (y_tilde(t) + theta1_vals * y_tilde(t-1) + delta_e) / u_tilde(t);
    theta2_min = (y_tilde(t) + theta1_vals * y_tilde(t-1) - delta_e) / u_tilde(t);
    
    plot(theta1_vals, theta2_max, 'r');
    plot(theta1_vals, theta2_min, 'b');

    colors = ['g', 'g', 'b', 'r'];
    fill([theta1_vals, fliplr(theta1_vals)], [theta2_max, fliplr(theta2_min)], colors(t), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

axis([-3 3 -3 3]);
xlabel('\theta_1');
ylabel('\theta_2');
title('Polytope Representation');
grid on, hold off
