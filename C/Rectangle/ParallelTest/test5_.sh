N=100
np=4
ol=10

printf "Fastest Newton LS on the subdomain: \n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'


printf "\nFAS on the subdomain:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'


printf "\nNGMRES with LS NPC on the subdomain:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This is fast. It's faster than the fastest newton linesearch.

# printf "\nQuasi-Newton on the subdomain:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_npc_snes_type fas -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## It's not fast.

printf "\nNewton TR on the subdomain:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type newtontr -sub_npc_snes_type fas -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This is pretty fast. It's on the same level as NGMRES.

printf "\n --- not using NASM beyond this point --- \n"


printf "\nNGMRES with Newton LS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -npc_snes_type newtonls -npc_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This is insanely fast.

printf "\nNGMRES with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -t1_problem ex10 -snes_type ngmres -npc_snes_type fas -npc_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
##
## This even faster.
