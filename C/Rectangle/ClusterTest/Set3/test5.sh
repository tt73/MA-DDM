#!/bin/bash

## These tests are for the NASM with a different nonlinear solver.

N=100
np=4
ol=10

alias format='-snes_converged_reason | grep "*Problem\|*Error\|WTime\|Nonlinear solve"'


printf "Newtonls with optimal settings: \n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1 format
## For comparison.

printf "\nTR:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type newtontr -sub_snes_trtol 1e-6 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type eisenstat format
## Converges without NPC
## Could be annoying to tune.

printf "\nFAS:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_fas_coarse_pc_type eisenstat format
## Converges without NPC

printf "\nAnderson with LS NPC:\n"
lsnpc=' -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1'
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc format

printf "\nAnderson with FAS NPC:\n"
fasnpc=' -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_pc_type eisenstat -sub_npc_fas_coarse_snes_max_it 1'
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type anderson -sub_snes_max_it 1 $fasnpc format

printf "\nNGMRES with LS NPC:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc format
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -sub_snes_max_it 1 -snes_converged_reason -snes_monitor -snes_view
## NGMRES doensn't converge without a preconditioner
## This is good

printf "\nNGMRES with FAS NPC:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_snes_max_it 1 $fasnpc format
## This is good

# printf "\nNGMRES with TR NPC:\n"
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_snes_max_it 1 -sub_npc_snes_type newtontr -sub_npc_pc_type eisenstat -sub_npc_snes_max_it 1 format
## This is bad. Trustregion is not a good preconditioner.


printf "\nQN with LS NPC:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 1 $lsnpc format
## It's not fast, but it has fewer DDM iterations.

printf "\nQN with FAS NPC:\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 1 $fasnpc format
## It's not fast, but it has fewer DDM iterations.

