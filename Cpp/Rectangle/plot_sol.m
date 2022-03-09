%% Notes

% Run this command: ./test1 -da_refine 2 -snes_monitor -snes_converged_reason -snes_type newtontr
% Then run this script

clear, close all

load_u % load petsc's numerical solution, u
load_exact % load exact solution coded in petsc
load_err

s = length(u); 
n = sqrt(s);

%% visualizing the solution
% The solution is one long M*N by 1 vector.
% We need to reshape it into a matrix before we plot it.

figure
u_grid = reshape(u,n,n);
surf(u_grid)
title('Numerical solution')

figure 
u_exac = reshape(u_exact,n,n);
surf(u_exac)
title('Exact solution')

figure
surf(u_grid-u_exac)
title('error')

fprintf('|error|_inf = %f',norm(u_grid(:)-u_exac(:),inf))

figure
err_s = reshape(err,n,n);
surf(err_s)