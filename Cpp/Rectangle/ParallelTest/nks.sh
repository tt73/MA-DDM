#!/bin/bash

# N = Number of interior points
N=20

# We test Newton-Krylov-Schwarz.

printf "\n\nNewton-Krylov-Schwarz\n"

printf "\n\nRun with 2 procs, default settings\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_view -snes_monitor -snes_type newtonls
# running this code with 2 procs automatically changes method to block-jacobi
# it works. we use newton's method with linesearch stepping and solve the Jacobian system with block-jacobi gmres

printf "\n\nConjugate Gradient\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_monitor -snes_type newtontr -snes_max_it 50 -ksp_type cg -snes_view
# this also works. however we have to change linesearch to trust region

printf "\n\nOther types of Krylov methods: \n"
../test1 -snes_type newtonls -help | grep ksp_type


# Trying out ksp types that do not give an error
printf "\n\nTrying out bunch of available ksp types\n"
ksps="pipecg pipecgrr pipelcg pipeprcg pipecg2 cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls"
for s in $ksps
do
   printf "\nTrying out $s:\n"
   mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type newtonls -ksp_type $s -log_view | grep 'error\|Nonlinear\|Time (sec):'
done
# they either diverge with 0 itersations or converges with 8-10 iterations 