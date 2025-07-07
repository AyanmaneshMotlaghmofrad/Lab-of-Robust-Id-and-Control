clc
clear variables
close all
format compact

%%
load data_exam_2A_hammer.mat

figure(1)
title('Input')
stem(Input,'LineWidth',1.5)
xlabel('sample numbers')
ylabel('The value of the input')
grid on,

figure(2)
title('Output')
stem(Output,'LineWidth',1.5)
xlabel('sample numbers')
ylabel('The value of the output')
grid on,
%% 
N = length(Input);
na = 1;
nb = 2;
tmin = max([na+1,nb]);
size_theta = 4;

y_tilde = Output;
u = Input;
delta_eta = 0.02;
%% Defining the equality constraints
for c = 1:N-tmin+1

    % support matrix  
    supp_mat = zeros(8,size_theta + N);
    
    %theta columns
    supp_mat(2:4,1:3) = eye(3);
    supp_mat(6:8,1:3) = eye(3);
    supp_mat(3:4,4) = ones(2,1);
    
    %eta columns
    supp_mat(5:6,size_theta+c:size_theta+c+1) = flip(eye(2));
    % coef vector
    coef_vec = [y_tilde(c+tmin-1) y_tilde(c+tmin-2) -u(c+tmin-1) -u(c+tmin-2)...
        -1 -1 -u(c+tmin-1)^3 -u(c+tmin-2)^3]';
    % equation structure
    ineqPolySys{c}.noTerms = 8;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.dimVar = size_theta + N;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vec;


end

%% POP parameters
ubd = [1e10*ones(4,1); delta_eta*ones(N,1)];
lbd = -ubd;

param.relaxOrder = 1;
param.POPsolver = 'active-set'

%% PUI lower bound

for i = 1:size_theta

    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.typeCone = 1;
    objPoly.dimVar = N + size_theta;

    supp = zeros(1, size_theta + N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = 1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_min(i) = POP.xVect(i);
    sol_refined_min(i) = POP.xVectL(i)
end

%%
% PUI upper bound
for i = 1:size_theta
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = size_theta + N;
    objPoly.typeCone = 1;

    supp = zeros(1, size_theta + N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i)
end

PUI = [sol_refined_min', sol_refined_max']
