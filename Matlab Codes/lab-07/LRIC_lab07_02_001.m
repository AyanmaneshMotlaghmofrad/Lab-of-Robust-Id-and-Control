clear variables
close all
clc
format compact

%% loading data

load data_lab_8.mat

y_tilde = ytilde;
clear ytilde
%% defining some constants
na = 1;
nb = 3;
tmin = max([na +1, nb]);

N = length (u);
size_theta = 3;
% assumption on the bound of the OE noise
delta_eta =0.01;

%% Defining the equality constraints

for c = 1: N-tmin+1
    
    % shaping the support matrix
    supp_mat = zeros(7,size_theta+N-1);

    % theta columns
    supp_mat(2:4,1:3) = eye(3);
    supp_mat(5:7,1:3) = [0 0 0
                         1 0 0
                         0 0 1];

    % eta columns
    supp_mat(5:7,c + size_theta:c + size_theta+1) = [0 1
                                                     1 0
                                                     1 0];
    
  
    % coefficient vector
    coef_vect = [y_tilde(c+tmin-1) -y_tilde(c+tmin-2) -u(c+tmin-2)^2 ...
        -y_tilde(c + tmin -2)*u(c + tmin - 3) -1  +1 u(c+tmin-3)]';
    
    ineqPolySys{c}.noTerms = 7;
    ineqPolySys{c}.dimVar = size_theta + N-1;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vect;

end


%% Defining SparsePOP parameters and bounds of the parameters.

ubd = [100*ones(3,1);delta_eta*ones(N-1,1)];
lbd = - ubd;

param.relaxOrder = 1;
param.POPsolver = 'interior-point';

%% SparsePOP for the lower bound

for i = 1:size_theta

    % support matrix
    supp = zeros(1,size_theta+N-1);
    supp(i) = 1;

    % 
    objPoly.noTerms = 1;
    objPoly.dimVar = size_theta + N-1;
    objPoly.degree = 1;
    objPoly.typeCone = 1;
    objPoly.supports = supp;
    objPoly.coef = 1;

    [a,b,POP] = sparsePOP(objPoly, ineqPolySys,lbd, ubd, param);
    sol_relaxed_min(i) = POP.xVect(i);
    sol_refined_min(i) = POP.xVectL(i);

end

%% SparsePOP solving for the maximum

for i = 1:size_theta

    % support matrix
    supp = zeros(1,size_theta+N-1);
    supp(i) = 1;

    % 
    objPoly.noTerms = 1;
    objPoly.dimVar = size_theta + N -1;
    objPoly.degree = 1;
    objPoly.typeConde = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [a,b,POP] = sparsePOP(objPoly, ineqPolySys,lbd, ubd, param);
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i);

end

PUI = [sol_refined_min' sol_refined_max']