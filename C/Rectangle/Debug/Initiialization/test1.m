% The goal of this test is to run the code for just 1 iteration of Newton
% and see how the initial guess evolves. 


N = 32; 
h = 2/(N+1);

cmd = sprintf("../../maddm -N %d -sol -debug -snes_max_it 1 -init_type coarse",N);

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

