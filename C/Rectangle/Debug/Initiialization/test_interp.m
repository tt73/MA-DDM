% The goal of this test is to run the code for just 1 iteration of Newton
% and see how the initial guess evolves. 


N = 16; 
k = 4;
N_cwb = round((N-1)/k) + 1

ratio = (N-1)/(N_cwb-1)



cmd = sprintf("../../interp -problem ex3 -N %d -coarseness %d",N,k);
out = system(cmd);

% load the petsc solutions 
load_cnb
load_cwb
load_interp

% reshape the solutions 
u_cnb     = reshape(u_cnb,sqrt(length(u_cnb)),sqrt(length(u_cnb)));
u_cwb     = reshape(u_cwb,sqrt(length(u_cwb)),sqrt(length(u_cwb)));
u_interp = reshape(u_interp,sqrt(length(u_interp)),sqrt(length(u_interp)));

figure 
subplot(131), surf(u_cnb), title('coarse solution without border')

subplot(132), surf(u_cwb), title('coarse solution with border') 

subplot(133), surf(u_interp), title('interpolated solution')


