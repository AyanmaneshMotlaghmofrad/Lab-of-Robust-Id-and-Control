% laboratory of robust identificatin
% lab 05_1_1
% Identification of two-input-one-output system with output-error noise
% structure
clc
clear variables
close all
format compact
rng('default')
%% Defining the system to be identified
% first transfer function
theta_1 = [1.4 0.76 0.184 2 1.8 0.68]
num1 = [0 theta_1(4:6)]
den1 = [1 theta_1(1:3)]
G1 = tf(num1,den1,-1,'var','z^-1')

%second transfer function
theta_2 = [1.4 0.76 0.184 3 3.4 1.28]
num2 = [0 theta_2(4:6)]
den2 = [1 theta_2(1:3)]
G2 = tf(num2,den2,-1,'var','z^-1')

%% simulation
N = 50; %simulation samples

u1 = 9*rand(N,1)-4;
u2 = 9*rand(N,1)-4;

% Noiseless output
w = lsim(G1,u1) + lsim(G2,u2);

% Noise assumption
Delta_eta = [-0.01 0.01];

% Noise corruption
eta = 2*Delta_eta(2)*rand(N,1)-Delta_eta(2);
y_tilde = w + eta;

%% Defining POP problem 

size_theta = 9; %considering a 2nd-order LTI system
na = 3;
nb = 3;
tmin = max([na+1,nb]);

%% main equation

for c = 1:N

      %initiation
    supp_mat = zeros(4,size_theta+3*N); 
% a1 a2 a3 b11 b21 b31 b12 b22 b32,[z1(1)..z1(N)],[z2(1)..z2(N)],[eta(1)..eta(N)]
    
    %rows of support matrix
    %y(1),-z1(1),-z2(1), -eta(1)
    
    % theta columns: not needed, all zero

    % z1 columns: (third row)
    supp_mat(2,c+size_theta) = 1;

    % z2 columns:
    supp_mat(3,c+size_theta+N) = 1;
    % eta columns:
    supp_mat(4,c+size_theta+2*N)=1;
    
    %coef vector
    coef_vec = [y_tilde(c) -1 -1 -1]';
    
    % equation structure
    ineqPolySys{c}.noTerms = 4;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.dimVar = size_theta + 3*N;
    ineqPolySys{c}.degree = 1;
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vec;
end

%% equality constraint z2
for c = N+1:2*N-tmin+1
    i = c-N;
    %initiation
    supp_mat = zeros(7,size_theta+3*N); 

% a1 a2 a3 b11 b21 b31 b12 b22 b32,[z1(1)..z1(N)],[z2(1)..z2(N)],[eta(4)..eta(N)]
    
    %rows of support matrix
%z2(k),+a1z2(k-1),a2z2(k-2),a3z2(k-3)],-b21u2(k-1),-b22u2(k-2),-b32(k-3)u2
    
% theta columns
    supp_mat(:,1:9) = [zeros(1,9);[[eye(3);zeros(3)],[zeros(6,3)],[zeros(3);eye(3)]]];
    
     % z1: not needed

     % z2: 
    supp_mat(1:4,i+N+size_theta:i+N+size_theta+tmin-1)=flip(eye(4));
  
    %coef_vec
    coef_vec = [1 1 1 1 -u2(i+tmin-2) -u2(i+tmin-3) -u2(i+tmin-4)]';
    % Simone: the correct index of u1 is as follows:
    %coef_vec = [1 1 1 1 -u2(i+tmin-1) -u2(i+tmin-2) -u2(i+tmin-3)]';

    ineqPolySys{c}.noTerms = 7;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.dimVar = size_theta + 3*N;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vec;

end

%% partial output 1 (z1)
for c = 2*N-tmin+2: 3*N-2*tmin+2
    i = c - (2*N-tmin+1); %dummy counter
    %initiation
    supp_mat = zeros(7,size_theta+3*N); 

% a1 a2 a3 b11 b21 b31 b12 b22 b32,[z1(1)..z1(N)],[z2(1)..z2(N)],[eta(4)..eta(N)]
    
    %rows of support matrix
%z1(k),+a1z1(k-1),a2z1(k-2),a3z1(k-3)], -b1u1(k-1),-b2u1(k-2),-b3(k-3)..
    % theta columns
    supp_mat(:,1:9) = [zeros(1,9);[eye(6),zeros(6,3)]];

    % z1
    supp_mat(1:4,i+size_theta:i+size_theta+tmin-1)=flip(eye(4));
    % z2: not needed!
    
    %coef vetor
    coef_vec = [1 1 1 1 -u1(i+tmin-2) -u1(i+tmin-3) -u1(i+tmin-4)]';
    % Simone: the correct index of u1 is as follows:
    %coef_vec = [1 1 1 1 -u1(i+tmin-1) -u1(i+tmin-2) -u1(i+tmin-3)]';
    
    % equation structure
    ineqPolySys{c}.noTerms = 7;
    ineqPolySys{c}.typeCone = -1; %equality
    ineqPolySys{c}.dimVar = size_theta + 3*N;
    ineqPolySys{c}.degree = 2;
    ineqPolySys{c}.supports = supp_mat;
    ineqPolySys{c}.coef = coef_vec;

end



%% 
% the bound of the noises can be consider as a lbd and ubd for the noise 
%variables

% lower bound
lbd = [-1e10*ones(size_theta+2*N,1);-Delta_eta(2)*ones(N,1)];
% upper bound
ubd = -lbd;

param.relaxOrder = 1;
param.POPsolver = 'active-set';


%% PUI lower bound
for i = 1:size_theta
    objPoly.noTerms = 1;
    objPoly.degree = 1;
    objPoly.dimVar = size_theta + 3*N;
    objPoly.typeCone = 1;

    supp = zeros(1, size_theta + 3*N);
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
    objPoly.dimVar = size_theta + 3*N;
    objPoly.typeCone = 1;

    supp = zeros(1, size_theta + 3*N);
    supp(i) = 1;
    objPoly.supports = supp;
    objPoly.coef = -1;

    [~,~,POP] = sparsePOP(objPoly, ineqPolySys, lbd, ubd, param);
    sol_relaxed_max(i) = POP.xVect(i);
    sol_refined_max(i) = POP.xVectL(i)
end

PUI = [sol_refined_min', sol_refined_max']