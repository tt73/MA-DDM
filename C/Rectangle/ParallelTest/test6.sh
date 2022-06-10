printf "\nNGMRES with Newton LS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -npc_snes_type newtonls -npc_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This is insanely fast.

printf "\nNGMRES with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -npc_snes_type fas -npc_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This even faster.
