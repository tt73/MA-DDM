#!/bin/bash

# N = Number of interior points
N=20

# We test Newton-Krylov-Schwarz.
#

printf "\n\nNewton-Krylov-Schwarz\n"

printf "\n\nRun with 2 procs, default settings\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_view -snes_monitor -snes_type newtonls
# running this code with 2 procs automatically changes method to block-jacobi
# it works. we use newton's method with linesearch stepping and solve the Jacobian system with block-jacobi gmres

printf "\n\nConjugate Gradient\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_monitor -snes_type newtontr -snes_max_it 50 -ksp_type cg -snes_view
# this also works. however we have to change linesearch to trust region