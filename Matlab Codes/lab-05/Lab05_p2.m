clear all
close all
clc

%% Problem Definition
z = tf('z',1);
theta_1 = [1.4, 0.76, 0.184, 2, 1.8, 0.68];
theta_2 = [1.4, 0.76, 0.184, 3, 3.4, 1.28];
G_p1 = tf(theta_1(4:6), [1, theta_2(1:3)], -1);
G_p2 = tf(theta_2(4:6), [1, theta_2(1:3)], -1);
G_p = [G_p1, G_p2];
size_theta = size(theta_1);

%% Simulation
N = 50;
u_1 = 9*rand([N,1]) -4;
u_2 = 9*rand([N,1]) -4;
u = [u_1, u_2];
w = lsim(G_p, u);


%% Data corruption 
Delta_eta = [-0.1, 0.1];
Delta_eps1 = [-0.001, 0.001];
Delta_eps2 = [-0.001, 0.001];
eta = 2*Delta_eta(2) .* rand(N,1) + Delta_eta(1);
eps1 = 2*Delta_eps1(2) .*rand(N,1) + Delta_eps1(1);
eps2 = 2*Delta_eps2(2) .*rand(N,1) + Delta_eps2(1);

y_tilde = w + eta;
u_tilde1 = u_1 + eps1;
u_tilde2 = u_2 + eps2;

%% PUI computation using sparsePOP
% equality constraint 1
for c = 1:N
    ineqPolySys{c}.noTerms = 4;
    ineqPolySys{c}.degree = 1;
    ineqPolySys{c}.dimVar = 2*size_theta(2)+5*N;
    ineqPolySys{c}.typeCone = -1; % equality

    % support matrix
    supp_eq1 = zeros(4, 2*size_theta(2)+5*N);
    supp_eq1(2, c+2*size_theta(2)+2*N) = 1;
    supp_eq1(3, c+2*size_theta(2)) = 1;
    supp_eq1(4, c+2*size_theta(2)+N) = 1;
    ineqPolySys{c}.supports = supp_eq1;

    % coefficients matrix
    coeffs_eq1 = [y_tilde(c), -1, -1, -1]';
    ineqPolySys{c}.coef = coeffs_eq1;
end

% equality constraint 2
for c = N+1:2*N-3
    i = c-N;

    ineqPolySys{c}.noTerms = 10;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.dimVar = 2*size_theta(2)+5*N;
    ineqPolySys{c}.typeCone = -1; % equality

    % support matrix
    supp_eq2 = zeros(10, 2*size_theta(2)+5*N);
    supp_eq2(2:7, 1:size_theta(2)) = eye(6);
    supp_eq2(1:4, i+2*size_theta(2):i+2*size_theta(2)+3) = flip(eye(4));
    supp_eq2(8:end, 4:6) = eye(3);
    supp_eq2(8:end, i+2*size_theta(2)+3*N:i+2*size_theta(2)+3*N+2) = flip(eye(3));
    ineqPolySys{c}.supports = supp_eq2;

    % coefficients matrix
    coeffs_eq2 = [1, 1, 1, 1, -u_tilde1(i+2), -u_tilde1(i+1), -u_tilde1(i),...
        1, 1, 1]';
    ineqPolySys{c}.coef = coeffs_eq2;
end

% equality constraint 3
for c = 2*N-2:3*N-6
    i = c-(2*N-3);

    ineqPolySys{c}.noTerms = 10;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.dimVar = 2*size_theta(2)+5*N;
    ineqPolySys{c}.typeCone = -1; % equality

    % support matrix
    supp_eq3 = zeros(10, 2*size_theta(2)+5*N);
    supp_eq3(2:7, size_theta(2)+1:2*size_theta(2)) = eye(6);
    supp_eq3(1:4, i+2*size_theta(2)+N:i+2*size_theta(2)+N+3) = flip(eye(4));
    supp_eq3(8:end, 10:12) = eye(3);
    supp_eq3(8:end, i+2*size_theta(2)+4*N:i+2*size_theta(2)+4*N+2) = flip(eye(3));
    ineqPolySys{c}.supports = supp_eq3;

    % coefficients matrix
    coeffs_eq3 = [1, 1, 1, 1, -u_tilde2(i+2), -u_tilde2(i+1), -u_tilde2(i),...
        1, 1, 1]';
    ineqPolySys{c}.coef = coeffs_eq3;
end

% lower bound
lbd = [-1e10*ones(2*size_theta(2)+2*N,1); -Delta_eta(2)*ones(N,1);...
    -Delta_eps1(2)*ones(N,1); -Delta_eps2(2)*ones(N,1)];
% upper bound
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';

% PUI lower bound
for i = 1:2*size_theta(2)
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = 2*size_theta(2) + 5*N;
    objPoly.typeCone = 1;

    supp = zeros(1, 2*size_theta(2) + 5*N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = 1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_min(i) = POP.xVect(i);
    sol_refined_min(i) = POP.xVectL(i);
end

% PUI lower bound
for i = 1:2*size_theta(2)
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = 2*size_theta(2) + 5*N;
    objPoly.typeCone = 1;

    supp = zeros(1, 2*size_theta(2) + 5*N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i);
end

PUI = [sol_refined_min', sol_refined_max']



