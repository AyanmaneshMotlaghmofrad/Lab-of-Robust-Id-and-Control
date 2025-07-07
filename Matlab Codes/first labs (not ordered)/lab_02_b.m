clear variable
close all
format compact
clc

%%
x1_lw = -15;
x1_up = 15;
x2_lw = -12;
x2_up = 12;

x1 = linspace(x1_lw,x1_up);
x2 = linspace(x2_lw, x2_up);

[x1,x2] = meshgrid(x1,x2);

z = x1.^2 + x2.^2 - 4.*x1 - 6.*x2 + 13;

%objective function
figure(1)
contour(x1,x2,z,5*(1:50)')
%%
%equality constraint

x2_1 = linspace(0.01,x2_up);
x2_2 = linspace(x2_lw,-0.01);
x1_1 = (10+7.*x2_1)./x2_1;
x1_2 = (10+7.*x2_2)./x2_2;
hold on;
plot(x1_1,x2_1,'k')
plot(x1_2,x2_2,'k')
axis([x1_lw,x1_up,x2_lw,x2_up])

% inequality constarints

x2_ineq = -0.5.*x1 +1;
plot(x1,x2_ineq,'b')

