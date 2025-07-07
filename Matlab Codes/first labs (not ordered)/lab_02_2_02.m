clc
clear variables
close all
format compact
rng('default')
%% defining the system
%LTI 2nd order, a-priori info
theta = [-0.7647,0.3012,0,32.24,21.41];

num = [theta(3) theta(4) theta(5)]
den = [1 theta(1) theta(2)];

Gp = tf(num,den,-1);

%% Simulation
N = 1000;
na = 2;
nb = 3;
tmin = max([na+1,nb]);

y = zeros(N,1);
u = rand(N,1);
delta_e = 0.01;
for t = tmin:N
    e = 2*delta_e*rand -delta_e; %uniform between[-delta_e,delta e]
    y(t) = -theta(1)*y(t-1) -theta(2)*y(t-2) + theta(3)*u(t) + theta(4)*u(t-1)+...
        theta(5)*u(t-2)+ e;
end

%% PUI problem solving by linprog
F = eye(5);
A1 = [-y(tmin-1:N-1), -y(tmin-2:N-2) u(tmin:N) u(tmin-1:N-1) u(tmin-2:N-2)];
A =[A1;-A1]; %for linprog
b1 = y(tmin:N) + delta_e;
b2 = -y(tmin:N) + delta_e;
b = [b1;b2];

% Uncomment if you want custom options for linprog:
% options = optimoptions('linprog', 'Display', 'none');

PUI = zeros(5,2);
sign = -1; % Initialization

for j = 1:2
    sign = -1*sign; % Alternate sign for each iteration
    for i = 1:5
        % Modify the linprog call to fit the correct syntax
        PUI(i,j) = F(i,:) * linprog(sign*F(i,:)', A, b);
    end
end
