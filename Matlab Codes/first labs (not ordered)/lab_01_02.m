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
N = 1000;
u = rand(N,1);
y = lsim(Gd,u); %y is the noise-free output of the system

%% Shaping regression matrix A using the simulation, 
% using the noise-free input and output samples.
% Assumption about the systerm
% LTI
% second-order => na = 2
na = 2;
nb = 3;
t_min = max(na+1, nb);
A = [-y(t_min-1:N-1) -y(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_LS = y(t_min:N);
% theta = [a1 a2 b1 b2 b3]
theta_LS = A \ y_LS
[num,den] = tfdata(Gd,'v')

theta_true = [den(2) den(3) num]';

%% Output error simulation
% or OE
e = 5*randn(N,1); %i.i.d noise
y_OE = lsim(Gd,u) + e;

%% Shaping regression matrix A using the simulation, 
% OE experiment
% Assumption about the systerm
% LTI
% second-order => na = 2
na = 2;
nb = 3;
t_min = max(na+1, nb)
A_OE = [-y_OE(t_min-1:N-1) -y_OE(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_LS_OE = y_OE(t_min:N);
% theta = [a1 a2 b1 b2 b3]
theta_OE = A_OE \ y_LS_OE
[num,den] = tfdata(Gd,'v')
%with 1000 samples it does not work


% Here the second assumption of consistency
% property of Least-Square does not hold.
%% Equation-error simulation

y_EE = lsim(Gd,u) + lsim(tf(1,den,1),e);

%% Shaping regression matrix A using the simulation, 
% EE experiment
% Assumption about the systerm
% LTI
% second-order => na = 2
na = 2;
nb = 3;

t_min = max(na+1, nb)
A_EE = [-y_EE(t_min-1:N-1) -y_EE(t_min-2:N-2) u(t_min:N) u(t_min-1:N-1) u(t_min-2:N-2)];
y_LS_EE = y_EE(t_min:N);
% theta = [a1 a2 b1 b2 b3]
theta_EE = A_EE \ y_LS_EE
[num,den] = tfdata(Gd,'v')
%with 1000 samples it does not work
%with 10000 it works, 

% Consistent with the second assumption of
% the consistency property

%% COMPARISON OF RESULT
Compare_THETA = [theta_true theta_EE theta_OE]
% Err_perc_EE = (theta_true-theta_EE)./theta_true*100;
% Err_perc_OE = (theta_true-theta_OE)./theta_true*100;

% Error_EE_and_OE = [Err_perc_EE Err_perc_OE]

%% EE norm infinity
A_EE_inf =[[-A_EE -one_vec];[A_EE -one_vec]];
b_EE_inf =[-y_EE(3:N);y_EE(3:N)];
F = eye(6);

PUI = zeros(5,2);
sign = -1; %initialization

for j = 1:2
    sign = -1*sign;
    for i = 1:5
        PUI(i,j) = F(i,:)*linprog(sign*[0,0,0,0,0,1]',A_EE_inf,b_EE_inf);
    end
end

%% norm infinity OE
one_vec = ones(N-na,1);
A_OE_inf =[[-A_OE -one_vec];[A_OE -one_vec]];
b_OE_inf =[-y_OE(3:N);y_OE(3:N)];
F = eye(6);

PUI = zeros(5,2);

PUI(:,1) = [1 1 1 1 1 0]*linprog(sign*[0,0,0,0,0,1]',A_OE_inf,b_OE_inf);
PUI(:,2) = [1 1 1 1 1 0]*linprog(sign*[0,0,0,0,0,-1]',A_OE_inf,b_OE_inf);
