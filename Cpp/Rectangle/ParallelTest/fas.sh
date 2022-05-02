#!/bin/bash

# N = Number of interior points
N=20

# We want to solve det(D^2u)=f using schwarz on the outside to decompose
# the domain and then use some non-linear solver on the local problems.
# The partitioning of the domain is automatic. Combining of the local
# solutions is also automatic. It is handled with snes type NASM.
printf "\n\nFull approximation scheme (FAS) options\n"
../test1 -snes_type fas -help | grep fas

printf "\n\nTest 1: Single proc \n"
mpiexec -n 1 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_view -snes_type fas
# it works. it does Newton linesearch with cubic backtracking. The linear solve is direct (LU decomp).

printf "\n\nTest 2: Two procs\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type fas -snes_max_it 50
# stagnation - dosen't converge to right solution
