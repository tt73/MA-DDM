% System Command

[status,cmd] = system('ls -l')

% MATLAB + PETSc
% Im reading this webpage: https://petsc4py.readthedocs.io/en/stable/manual/matlab/



% status = system('mpiexec ./poisson -mat_view :A.m:ascii_matlab -vec_view :b.m:ascii_matlab')

[status,cmd] = system('which mpiexec')


[~,pth] = system('echo $PATH')

pth = strsplit(pth,':')



[status,cmd] = system('mpiexec ./poisson -h')

[~,pth] = system('echo $PATH')

pth = strsplit(pth,':')