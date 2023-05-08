% The goal of this test is to run the code for just 1 iteration of Newton
% and see how the initial guess evolves. 


N = 17; 
k = 2;

cmd = sprintf("../../maddm -problem ex3 -N %d -init_type coarse -coarseness %d -sol -debug",N,k);
out = system(cmd);


% load the petsc solutions 
load_u
load_initial
load_interp

% reshape the solutions 
u     = reshape(u,sqrt(length(u)),sqrt(length(u)));
u_initial     = reshape(u_initial,sqrt(length(u_initial)),sqrt(length(u_initial)));
u_interp = reshape(u_interp,sqrt(length(u_interp)),sqrt(length(u_interp)));

figure 
subplot(131), surf(u), title('solution')

subplot(132), surf(u_initial), title('initial guess') 

subplot(133), surf(u_interp), title('interpolated solution')


