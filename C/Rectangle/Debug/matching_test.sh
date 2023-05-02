## Goal:
# Want to compare the
#
#


# load petsc
module load gnu8 mpich petsc/3.12.0

clear

## Settings
N=64     ## N=2^5
prob=ex1 ## gaussian [-1,1]^2
h=$(echo "scale = 8; 2/($N+1)" | bc)
h2=$(echo "scale = 8; $h*$h" | bc)
np=4     ## four subdomains

## Newton
# 1 iteration only
# basic linesearch

## linear solver
# use LU factorization to solve directly
# ksp = preonly
# pc = lu

## Overlap
# By percentage 10% overlap
# When np=4 and N=32, this is just ceil(0.1*32/2) = 2.

## run with rtol
# printf "Running with relative residue tolerance\n"
# mpiexec -np 4 ../maddm -N $N -problem $prob -sin -snes_converged_reason -snes_rtol 1e-8 -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic #-snes_linesearch_damping 0.5 -sub_snes_monitor

## run with atol to h^2
# printf "\n\n\nRunning with absolute residue tolerance to h^2 = $h2\n"
# mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h2 -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic #-snes_linesearch_damping 0.5


## SIN test...
## Does not match with matlab code
##
## run with atol to h
# printf "\n\n\nRunning with absolute residue tolerance to h = $h\n"
# printf "mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic  -snes_nasm_damping 0.9\n"
# mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic -snes_nasm_damping 0.9


## run with atol to h (BT on)
printf "\n\n\nRunning with absolute residue tolerance to h = $h\n"
printf "mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type bt  -snes_nasm_damping 0.9\n"
mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type bt





## one more run with N = 128
# N=128
# h=$(echo "scale = 8; 2/($N+1)" | bc)
# printf "\n\n\nRunning with absolute residue tolerance to h = $h\n"
# mpiexec -np 4 ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic -sub_snes_monitor



## run with atol to h for "exact Newton"
# printf "\n\n\nRunning with absolute residue tolerance to h = $h\n"
# printf "mpiexec -np 4 ../maddm -N $N -problem ex1 -htn -snes_converged_reason -snes_rtol 1e-99 -snes_atol $h -sub_snes_atol $h -sub_snes_rtol 1e-99 -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic -snes_max_it 100
# \n"
# mpiexec -np 4 ../maddm -N $N -problem ex1 -htn -snes_converged_reason -snes_rtol 1e-99 -snes_atol $h -sub_snes_force_iteration -sub_snes_atol $h -sub_snes_rtol 1e-99 -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type bt -snes_max_it 100



## run with atol to h for "exact Newton"
# printf "\n\n\nRunning with absolute residue tolerance to h = $h\n"
# printf "mpiexec -np 4 ../maddm -N $N -problem ex1 -htn -snes_converged_reason -snes_rtol 1e-99 -snes_atol $h -sub_snes_atol $h -sub_snes_rtol 1e-99 -snes_monitor -snes_view -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic -snes_max_it 100
# \n"
# mpiexec -np 4 ../maddm -N $N -problem ex1 -htn -snes_converged_reason -snes_rtol 1e-99 -snes_atol $h -sub_snes_force_iteration -sub_snes_atol $h -sub_snes_rtol 1e-99 -snes_monitor -snes_view -sub_ksp_type gmres -sub_snes_linesearch_type bt -snes_max_it 100
