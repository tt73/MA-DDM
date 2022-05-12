%% Notes

% Run this command: ./test1 -da_refine 5 -t1_dim 1
% To run it in parallel 
% Then run this script

clear, close all

load_u % load petsc's numerical solution, u
load_exact % load exact solution coded in petsc

n = length(u);


%% visualizing the solution

% assume x = 0:h:1
h = 1/(n-1);
x = 0:h:1;

figure
plot(x,u,'o-')
hold on 
plot(x,u_exact,'x-');
legend('Numerical solution','Exact solution')

figure
plot(x,u-u_exact)
title('error')

fprintf('|error|_inf = %f\n',norm(u-u_exact,inf))
