% System Command

status = system('ls -l')

% MATLAB + PETSc
% Im reading this webpage: https://petsc4py.readthedocs.io/en/stable/manual/matlab/

status = system('mpiexec ./poisson -mat_view :A.m:ascii_matlab -vec_view :b.m:ascii_matlab')

