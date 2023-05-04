% The goal of this test is to run the code for just 1 iteration of Newton
% and see how the initial guess evolves. 


N = 31; 
h = 2/(N+1);

init = "coarse";
% init = "zeros";
% init = "random";
cmd = sprintf("../../maddm -problem ex2 -N %d -sol -debug -snes_max_it 100 -init_type %s -snes_atol 1e-10",N,init);

out = system(cmd);

load_initial % run petsc-generated code - creates `u_initial` variable 

load_u % run petsc-generated code - creates `u` variable

% reshape N^2 by 1 to N by N 
u_initial = reshape(u_initial,N,N);
u = reshape(u,N,N);


figure 
subplot(121) 
surf(u_initial) 
subplot(122)
surf(u)


N4h =@(N) (N+1)/4 - 1;
Nh =@(N4h) 4*(N4h+1)-1;

