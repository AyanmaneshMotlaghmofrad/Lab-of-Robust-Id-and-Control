clear variables
close all
clc
format compact
rng('default')
%% system description

s = tf('s');
Gp = 100/(s^2 + 1.2*s +1);

%% discretizing

Gd = c2d(Gp,1,'zoh');

%% tfdata

[z,p,k] = zpkdata(Gd,'v');

%% Noise-free simulation
N = 10000;
u = rand(N,1);
y = lsim(Gd,u); % w is the output of the system

%% Shaping regression matrix A using the simulation, 
% using the noise-free input and output samples.
% Assumption about the systerm
% LTI
% second-order => na = 2
eps = 5; %bound of error
%norm infinity
na = 2;
nb = 3;
t_min = max(na+1, nb)
A = [-y(t_min-1:N-1) -y(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_inf = y(t_min:N);
% theta = [a1 a2 b1 b2 b3]
for i = 1:5
    theta_inf(i) = A(:,i)'*(y_inf - eps*sign(A(:,i)))
end

[num,den] = tfdata(Gd,'v')

%% Error in variable simulation
% or EIV
e = 5*randn(N,1); %i.i.d noise
y_EIV = lsim(Gd,u) + e;

%% Shaping regression matrix A using the simulation, 
% EIV experiment
% Assumption about the systerm
% LTI
% second-order => na = 2
na = 2;
nb = 3;
t_min = max(na+1, nb)
A_EIV = [-y_EIV(t_min-1:N-1) -y_EIV(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_inf_EIV = y_EIV(t_min:N);
% theta = [a1 a2 b1 b2 b3]
for i = 1:5
    theta_inf_EIV(i) = A_EIV(:,i)'*(y_inf_EIV - eps*sign(A_EIV(:,i)))
end
theta_inf_EIV
[num,den] = tfdata(Gd,'v')
%with 1000 samples it does not work
%with 10000 it works, 

% Consistent with the consistency property
% we expected this result, since it satisfies 
% both assumptions of the consistency property.
%% Output-error simulation

y_OE = lsim(Gd,u) + lsim(tf(1,den,-1),e);

%% Shaping regression matrix A using the simulation, 
% OE experiment
% Assumption about the systerm
% LTI
% second-order => na = 2
na = 2;
nb = 3;
t_min = max(na+1, nb)
A_OE = [-y_OE(t_min-1:N-1) -y_OE(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_inf_OE = y_OE(t_min:N);
% theta = [a1 a2 b1 b2 b3]
for i = 1:5
    theta_inf_OE(i) = A_OE(:,i)'*(y_inf_OE - eps*sign(A_OE(:,i)))
end

[num,den] = tfdata(Gd,'v')
%with 1000 samples it does not work
%with 10000 it works, 

% Consistent with the consistency property
% we expected this result, since it satisfies 
% both assumptions of the consistency property.

%%