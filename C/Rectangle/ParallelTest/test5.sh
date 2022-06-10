## These tests are for the NASM with a different nonlinear solver.

N=100
np=4
ol=10

## The variable below contains commands to print the code output neatly.
## At the end of the run command, add the line "$p".
p="-snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'"


printf "Newtonls with optimal settings: \n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1 "$p"
## For comparison.

printf "\nFAS:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_fas_coarse_pc_type eisenstat "$p"
## FAS requires you to choose a SNES on the coarse level.
## Its not that fast.


printf "\nFAS with LS NPC:\n"
npcls='-sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1'
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_fas_coarse_pc_type eisenstat -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 "$p"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_fas_coarse_pc_type eisenstat -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -snes_view
## Even with a LS NPC it doesn't run that fast.

printf "\nNGMRES with LS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -sub_snes_max_it 1 "$p"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -sub_snes_max_it 1 -snes_converged_reason -snes_monitor -snes_view
## NGMRES doensn't converge without a preconditioner
## Faster than LS!

printf "\nNGMRES with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_pc_type eisenstat -sub_npc_fas_coarse_snes_max_it 1 -sub_snes_max_it 1 "$p"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type ngmres -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_pc_type eisenstat -sub_npc_fas_coarse_snes_max_it 1 -sub_snes_max_it 1 -snes_converged_reason -snes_view
#
#

# printf "\nQN:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type l2 -sub_npc_snes_type fas -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type l2 -snes_converged_reason -snes_view -snes_monitor
## QN dosen't converge without a NPC. It doesn't matter what linesearch is used in the QN.
## QN dosen't seem to converge with a LS NPC either.

printf "\nQN with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type basic -sub_snes_linesearch_max_it 1 -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_snes_max_it 1 -snes_converged_reason "$p"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type qn -sub_snes_linesearch_type basic -sub_snes_linesearch_max_it 1 -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_snes_max_it 1 -snes_converged_reason -snes_view -snes_monitor
## It's not fast, but it has fewer DDM iterations.


printf "\nTR with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type newtontr -sub_npc_snes_type fas -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_snes_type newtontr -sub_npc_snes_type fas -sub_npc_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_type gmres -sub_pc_type eisenstat -snes_converged_reason -snes_view -snes_monitor
## This is pretty fast. It's on the same level as NGMRES.
