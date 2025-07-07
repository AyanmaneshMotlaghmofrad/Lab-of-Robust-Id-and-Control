clc
clear variable 
close all
format compact
rng('default')

%% System definition
theta = [-0.7647,0.3012,0,32.24,21.41];
par_lenth = length(theta);
num = [theta(3) theta(4) theta(5)]
den = [1 theta(1) theta(2)];

Gp = tf(num,den,-1,'var','z^-1')

%% Simulation
N = 200; % to be varied
na = 2;
nb = 3;
tmin = max([na+1,nb]);

y = zeros(N,1);
u = rand(N,1);

%noise free simulation
for t = tmin:N
    y(t) = -theta(1)*y(t-1) -theta(2)*y(t-2) + theta(3)*u(t) + theta(4)*u(t-1)+...
        theta(5)*u(t-2);
end

%% Corrupting the data
%noise assumptions
delta_xi = 0.01; %to be varied
delta_eta = 0.01; % to be varied

%corrupting the input and output samples
xi = 2*delta_xi*rand(N,1) - delta_xi; 
u_tilde = u + xi;

eta = 2*delta_eta*rand(N,1) - delta_eta;
y_tilde = y + eta;


%% the equality constraint matrices
%shaping the support matrix
theta_sup = zeros(12,par_lenth);
theta_block = [1;1];

j = 1;
for i = 3:2:11
   theta_sup(i:i+1,j)= theta_block;
   j = j + 1 ;
end


%% eta rows (1-6)
eta_block=[0  0  0;
           0  0  1;
           0  0  0;
           0  1  0;
           0  0  0;
           1  0  0];
eta_sup= zeros(6,2*N,N-tmin);
eta_sup(:,:,1) = [eta_block, zeros(6,2*N-tmin)];
for i = 2:N-tmin+1
    eta_sup(:,:,i) = [zeros(6,i-1),eta_block,zeros(6,2*N-tmin-(i-1))];
end

%% xi rows (7-12)
xi_block=[0  0  0;
           0  0  1;
           0  0  0;
           0  1  0;
           0  0  0;
           1  0  0];
xi_sup= zeros(6,2*N,N-tmin);
xi_sup(:,:,N-tmin) =[zeros(6,2*N-tmin),xi_block];
for i = 1:N-tmin+1
    xi_sup(:,:,i) = [zeros(6,i-1),xi_block,zeros(6,2*N-tmin-(i-1))];
end

%% support matrix for equality
sup_mat = zeros(12,par_lenth + 2*N,N-tmin+1);
for i = 1:N-tmin+1
    sup_mat(:,:,i) = [theta_sup,[eta_sup(:,:,i);xi_sup(:,:,i)]];
end
%%
for i = 1:N-tmin+1
    %inequality constraints
    eq_coef = [y_tilde(tmin+i-1), -1,y_tilde(tmin+i-2),-1,y_tilde(tmin+i-3),-1,-u_tilde(tmin+i-1),+1,-u_tilde(tmin+i-2),1,-u_tilde(tmin+i-3),1]';
    ineqPolySys{i}.typeCone = -1;
    ineqPolySys{i}.dimVar = par_lenth + 2*N;
    ineqPolySys{i}.degree = 2;
    ineqPolySys{i}.noTerms = 12;
    ineqPolySys{i}.supports = sup_mat(:,:,i);
    ineqPolySys{i}.coef = eq_coef;
end
%% solving the optimization problem

%lower and upper bounds
lbd = [-1e10*ones(5,1);-delta_eta*ones(N,1);-delta_xi*ones(N,1)];
    
ubd = -lbd;

%relaxation Order
param.relaxOrder = 2;
param.POPsolver = 'interior_point';

%defining objPoly

obj_mat = [eye(par_lenth),zeros(par_lenth,2*N);eye(par_lenth),zeros(par_lenth,2*N)];
obj_coef = [ones(par_lenth,1);-ones(par_lenth,1)];

%% SparsePOP
for i = 1:2*par_lenth
    %the main equation
    objPoly{i}.typeCone = 1;
    objPoly{i}.dimVar = par_lenth+2*N;
    objPoly{i}.noTerms = 1;
    objPoly{i}.degree = 1;
    objPoly{i}.supports = obj_mat(i,:);
    objPoly{i}.coef = obj_coef(i);
    
    [a,b,POP] = sparsePOP(objPoly{i},ineqPolySys,lbd,ubd,param);
    sol_relaxed(i) = POP.xVect;
    sol_refined(i) = POP.xVectL;
end

%% comparing the results

