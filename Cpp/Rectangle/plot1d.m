%% Notes

% Run this command: ./test1 -da_refine 5 -t1_dim 1
% To run it in parallel 
% Then run this script

clear, close all

load_u % load petsc's numerical solution, u
load_exact % load exact solution coded in petsc

n = length(u); 

%% visualizing the solution
% The solution is one long M*N by 1 vector.
% We need to reshape it into a matrix before we plot it.

figure
plot(u)
title('Numerical solution')

figure 
plot(u_exact);
title('Exact solution')

figure
surf(u-u_exact)
title('error')

fprintf('|error|_inf = %f\n',norm(u-u_exact,inf))
