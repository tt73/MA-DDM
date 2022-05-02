#!/bin/bash

# N = Number of interior points
N=20

# We want to solve det(D^2u)=f using schwarz on the outside to decompose
# the domain and then use some non-linear solver on the local problems.
# The partitioning of the domain is automatic. Combining of the local
# solutions is also automatic. It is handled with snes type NASM.
printf "\n\nNon-linear additive schwarz options\n"
../test1 -help | grep nasm

printf "\n\nTest 1: One domain\n"
mpiexec -n 1 ../test1 -t1_N ${1:-$N} -snes_view -snes_type nasm
# it works. it does Newton linesearch with cubic backtracking. The linear solve is direct (LU decomp).

printf "\n\nTest 2: Two domains\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -snes_max_it 50
# stagnation - dosen't converge to right solution

printf "\n\nTest 3: Two domains, custom local domain solving\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -snes_max_it 50 -sub_snes_linesearch_type basic -sub_ksp_type gmres -sub_pc_type ilu
# stagnation again

printf "\n\nTest 4: Try -npc nasm \n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type newtonls -npc_snes_type nasm -snes_max_it 50
# still doesn't work

printf "\n\nTest 5: Try -npc nasm with simple stepping \n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -npc_snes_type nasm -snes_max_it 50 -npc_sub_snes_linesearch_type basic -snes_view
# still doesn't work
# NPC is very confusing, little documentation

printf "\n\nTest 6: Try the final jacobian setting\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -snes_max_it 50 -snes_nasm_finaljacobian_type finalinner -snes_view -snes_nasm_finaljacobian
# no

