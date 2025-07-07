close all
clear variables
clc
format compact

%% Defining system
z = tf('z',1);
theta = [-0.5 2];
Gp = theta(2)/(1 + theta(1)*z)

u_tilde = [4 -3 2 1]';
e = [0.05 -0.25 0.3 -0.5]';

%% simulation
% y_tilde(k) = -theta_1*y_tilde(k-1) + theta_2*u_tilde(k) + e(k)
na = 1; %the number of state variables or last samples of y
nb = 1;%the number of past inputs
N = 4;
t_min = max([na+1 , nb]);

y_tilde= zeros(4,1);
%initial value

for k = t_min:N
    y_tilde(k) = -theta(1)*y_tilde(k-1) + theta(2)*u_tilde(k) + e(k);
end

%% upper and lower band 
delta_e = 0.5;
y_up = y_tilde + 0.5;
y_low = y_tilde - 0.5;

plot(1:N, [y_up, y_tilde,y_low])

%% 
theta_1 = linspace(-5,5,1000)
figure,
hold on
for k = t_min : N
    k
    theta_2_up = (+delta_e + y_tilde(k) +theta_1*y_tilde(k-1))/u_tilde(k);
    theta_2_low =(-delta_e + y_tilde(k) + theta_1*y_tilde(k-1))/u_tilde(k);
    plot(theta_1,theta_2_low,'b')
    plot(theta_1,theta_2_up,'r')
end



%% stupid way of doing it
figure,
hold on
for theta_1 = -100:0.01:100
    for theta_2 = -100:0.01:100
        for k = t_min:N
            f = 1; %flag
            if(abs(y_tilde(k)+theta_1*y_tilde(k-1)-theta_2*u_tilde(k))>delta_e)
                f = 0;
                break
            end
        end
        if (f == 1)
            plot(theta_1,theta_2,'+');
        end
    end
end
