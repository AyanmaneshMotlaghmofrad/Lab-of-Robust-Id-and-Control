clear all
close all
clc

%% Problem Definition
z = tf('z',1);
theta = [-0.7647, 0.3012, 0, 32.24, 21.41];
G_p = tf(theta(3:5), [1, theta(1:2)], 1);
size_theta = size(theta);

%% Simulation
N = 100;
u = rand(N,1);
w = lsim(G_p, u);

%% Data corruption 
Delta_eta = [-1, 1];
Delta_eps = [-0.01, 0.01];

eta = 2*Delta_eta(2) .* rand(N,1) + Delta_eta(1);
eps = 2*Delta_eps(2) .* rand(N,1) + Delta_eps(1);

y_tilde = lsim(G_p, u) + eta;
u_tilde = u + eps;

%% PUI computation using sparsePOP
% equality constraint
for c = 1:N-2

    % equality support matrix
    %initialization
    supp_eq = zeros(12, size_theta(2)+2*N);

    %theta rows
    supp_eq(2:6, 1:size_theta(2)) = eye(size_theta(2)); 
    supp_eq(8:end, 1:size_theta(2)) = eye(size_theta(2));
    
    % eta rows
    supp_eq(7:9, size_theta(2)+c:size_theta(2)+c+2) = flip(eye(3));

    % xi rows
    supp_eq(10:12, size_theta(2)+c+N:size_theta(2)+c+2+N) = flip(eye(3));
    
    %equality coefficient vector
    coeffs_eq = [y_tilde(c+2), y_tilde(c+1), y_tilde(c),...
        -u_tilde(c+2), -u_tilde(c+1), -u_tilde(c), -1, -1, -1, 1, 1, 1]';
    
    % creating the structures
    ineqPolySys{c}.noTerms = 12;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.dimVar = size_theta(2) + 2*N;
    ineqPolySys{c}.typeCone = -1; % equality
    ineqPolySys{c}.supports = supp_eq;
    ineqPolySys{c}.coef = coeffs_eq;
end

% the bound of the noises can be consider as a lbd and ubd for the noise 
%variables

%% lower bound
lbd = [-1e10*ones(size_theta(2),1); -Delta_eta(2)*ones(N,1); -Delta_eps(2)*ones(N,1)];
% upper bound
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

%% PUI lower bound
for i = 1:size_theta(2)
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = size_theta(2) + 2*N;
    objPoly.typeCone = 1;

    supp = zeros(1, size_theta(2) + 2*N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = 1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_min(i) = POP.xVect(i);
    sol_refined_min(i) = POP.xVectL(i);
end

% PUI upper bound
for i = 1:size_theta(2)
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = size_theta(2) + 2*N;
    objPoly.typeCone = 1;

    supp = zeros(1, size_theta(2) + 2*N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i);
end

PUI = [sol_refined_min', sol_refined_max']
