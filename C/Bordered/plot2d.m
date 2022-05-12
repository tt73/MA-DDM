%% Notes
% To set dimension (1,2,3), run ./test1 -t1_dim 2
% But 2D is default so you don't really need to use the flag and 3D does
% not work atm. 

% To choose problem, run ./test1 -t1_problem ex10
% ex11 and ex12 are also available.

% For parallel, run mpiexec -n 4 ./test1

% To use a Jacobian free method, run ./test1 -snes_mf or -snes_fd

% Not all Nonlinear solvers converge. The ones that converge are aspin and
% newtontr. Run ./test1 -snes_type aspin or -snes_type newtontr

clear, close all

load_u % load petsc's numerical solution, u
load_exact % load exact solution coded in petsc

s = length(u);
n = sqrt(s);

%% Mesh 
% domain [0 1] x [0 1] 
% h = 1/(n-1), endpoints included

h = 1/(n-1); 
x = 0:h:1;
y = x; 
[X,Y] = meshgrid(x,y); 

%% Exact Solution 

% choose F
choice = 1;
switch(choice)
   case 1
      U = @(x,y) (exp((x.^2+y.^2)/2));
   case 2
      pos = @(x) max(x,0);
      U = @(x,y) (.5*pos(((x-.5).^2+(y-.5).^2).^.5-.2).^2);
   case 3
      U = @(x,y) (-sqrt(2-(x.^2+y.^2)));
end

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

fprintf('|error|_inf = %f\n',norm(u_grid(:)-u_exac(:),inf))

%% Compare 
figure 
title('Compare PETSc sol with MATLAB') 
h = surf(u_grid);
hold on 
e = surf(U(X,Y))