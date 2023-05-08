% The goal of this test is to run the code for just 1 iteration of Newton
% and see how the initial guess evolves. 


N = 17; 
k = 2;

cmd = sprintf("../../maddm -problem ex1 -N %d -init_type coarse -coarseness %d -sol -debug",N,k);
out = system(cmd);


% load the petsc solutions 
load_u
load_initial
load_interp

% reshape the solutions 
u         = reshape(u,sqrt(length(u)),sqrt(length(u)));
u_initial = reshape(u_initial,sqrt(length(u_initial)),sqrt(length(u_initial)));

figure 
subplot(121), surf(u), title('solution')

subplot(122), surf(u_initial), title('initial guess') 



