%% Notes
%
% I want to run petsc code and use Matlab for debugging and plotting.
%
% Im reading this webpage: https://petsc4py.readthedocs.io/en/stable/manual/matlab/

clear


%% System Command
% You can send commands to the shell with the `system` function.
% This line of code
[code,out] = system('ls -l');
disp("Result from ls -l command:")
disp(out)


%% Issue running programs
% For some reason, calling the program with the system function doesn't work.
%
M = 100;
N = 40;
cmd = sprintf('./poisson -da_grid_x %d -da_grid_y %d',M,N);
[code,out] = system(cmd);
% I get an error: './poisson: /usr/local/MATLAB/R2020b/sys/os/glnxa64/libgfortran.so.5: version `GFORTRAN_9' not found (required by /lib/x86_64-linux-gnu/libmpichfort.so.12)

%% The simplest work-around
% Quickest way of interfacing with Matlab is to run the code outside of Matlab.
% Then you can run this script to load the system.
%
% The poisson program is configured to print out the system and the
% solution in matlad syntax. The three scripts: load_mat.m, load_u.m, and
% load_b.m are generated automatically with PETSc.

% Run the program with ./poisson -da_grid_x M -da_grid_y N
% where M and N are integers of your choosing
M = 100
N = 50

% The following lines loads variables used in the petsc program

load_mat % load the A matrix

load_u % load petsc's numerical solution, u

load_b % load the rhs b vector

%% visualizing the solution
% The solution is one long M*N by 1 vector.
% We need to reshape it into a matrix before we plot it.

figure
spy(A)

figure
u_grid = reshape(u,M,N);
surf(u_grid)

