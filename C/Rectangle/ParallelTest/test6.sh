#!/bin/bash

N=100
np=4
ol=10

alias format='-snes_converged_reason | grep "*Problem\|*Error\|WTime\|Nonlinear solve"'

printf "\nLS:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2 format

printf "\nFAS:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type fas -fas_coarse_snes_linesearch_type l2 format

# printf "\nNGMRES:\n"
# mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -snes_ngmres_select_type difference


printf "\nNGMRES with Newton LS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_type basic -npc_snes_max_it 1 -snes_view


printf "\nNGMRES with FAS NPC:\n"
fasnpc='-npc_snes_type fas -npc_fas_coarse_snes_linesearch_type l2 -npc_snes_max_it 1 -npc_fas_levels_snes_type newtonls -npc_fas_levels_snes_max_it 6'
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres $fasnpc -snes_view


printf "\nFAS * NGMRES composite:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type composite -snes_composite_type multiplicative -snes_composite_sneses fas,ngmres  -snes_view
