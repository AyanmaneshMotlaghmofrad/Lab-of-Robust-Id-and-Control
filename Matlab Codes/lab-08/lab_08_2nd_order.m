% Third order controller 
% it seems like it does not guarantee the stability of the system
clc
clear variables 
close all
format compact


%% System and reference system
s = tf('s');
Gp = 4500/(s^2 + 16*s + 4500);

M = (51.54*s+861.8)/(s^2 + 56.97*s + 861.8);
%% simulation for obtaining Ts
t = linspace(1,2,10000);
u = ones(1,length(t));
y = lsim(Gp,u,t);
y_r = lsim(M,u,t);
figure(1)
title('step response')
plot(t,y,'LineWidth',1.5)
hold on
plot(t,y_r,'-r','LineWidth',1.5)

% settelling time = 1.5
N_trans_sampling = 75; % # of samples to capture the transient 50
Ts = 1.5/N_trans_sampling


%% discretizing the system
Ts = 0.02; %rounding

Gp_d = c2d(Gp,Ts,'zoh')
M_d = c2d(M,Ts,'zoh')

%% simulating experimental data

N = 50; % the number of samples

u = 3*rand(N,1)-1; %input

y = lsim(Gp_d,u); %noise-free data

delta_eta = 0.001;
eta = 2*delta_eta*rand(N,1)-delta_eta;

y_tilde = y + eta;

%% simulating S(t)

s = lsim(minreal(M_d/(1-M_d),1e-4),u);

%% SM problem for a third order controller K third order
na = 2;
nb = 2;
nk = 0;
tmin = max([na+1,nb+nk]);
size_rho = na + nb + 1;

for c = 1:N-tmin +1
    
    % support matrix
    supp_mat = zeros(3*(na+1),size_rho+N);

    %rho columns
    supp_mat(2:6,1:5) = eye(5);
    supp_mat(7:9,3:5) = eye(3);

    %eta columns
    supp_mat(7:9,size_rho + c:size_rho + c + 2) = flip(eye(3));

    % coef vector
    coef_vec = [s(c+tmin-1) s(c+tmin-2) s(c+tmin-3) ...
        y_tilde(c+tmin-1) y_tilde(c+tmin-2) y_tilde(c+tmin-3)...
        1 1 1]';
    % defining the equality constraints
    ineqPolySys{c}.noTerms = 3*(na+1);
    ineqPolySys{c}.dimVar = size_rho + N;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vec;

end

%% optimization parameters and bounds on the variables

ubd = [1e10*ones(size_rho,1);delta_eta*ones(N,1)];
lbd = -ubd;

param.relaxOrder = 1;
param.POPsolver = 'interior-point';

%% solving the optimization problem using SparsePOP
% lower bound
for i = 1:size_rho
    %support matrix
    supp = zeros(1,size_rho+N);
    supp(i) = 1;

    objPoly.noTerms = 1;
    objPoly.dimVar = size_rho + N;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = supp;
    objPoly.coef = 1;

    [a,b,POP] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    sol_relaxed_min(i) = POP.xVect(i);
    sol_refined_min(i) = POP.xVectL(i);

end

% upper bound

for i = 1:size_rho
    %support matrix
    supp = zeros(1,size_rho+N);
    supp(i) = 1;

    objPoly.noTerms = 1;
    objPoly.dimVar = size_rho + N;
    objPoly.typeCone = 1;
    objPoly.degree = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [a,b,POP] = sparsePOP(objPoly,ineqPolySys,lbd,ubd,param)
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i);

end

%%
CPUI = [sol_refined_min', sol_refined_max']

%% shaping the controller
rho_c = mean(CPUI,2);

Kc = tf([rho_c(3:end)'], [1 rho_c(1:2)'],Ts)

%% simulation
T_c = minreal((Kc*Gp_d)/(1+Kc*Gp_d),1e-4);%complementary sensitivity function

figure(1)
title('step response')
step(M_d)
hold on,
step(T_c)
legend('reference model','model with the designed controller')