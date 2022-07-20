#!/bin/bash -l
## We still use NASM as the gloabal method.
## We use a local solver other than Newton's method.

## NAMING
#SBATCH -J j5
#SBATCH -p public
#SBATCH -o slurmout5
#SBATCH -e slurmout5

## partition/queue
#SBATCH -p dms-cpu

## EMAIL NOTIFICATION
#SBATCH --mail-user tt73@njit.edu
#SBATCH --mail-type=END
#SBATCH -A tt73

## RUNTIME HOURS:MIN:SEC and MEMORY
#SBATCH -t 8:0:0
#SBATCH --mem=16G
#SBATCH -N 4

module load gnu8 mpich petsc

rm -f out5

N=200
np=4

printf "SIN:\n" >> out5
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out5

printf "TR:\n" >> out5
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type newtontr -sub_snes_trtol 1e-6 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type newtontr -sub_snes_trtol 1e-6 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type newtontr -sub_snes_trtol 1e-6 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_type newtontr -sub_snes_trtol 1e-6 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out5
## Its slow but it converges without NPC
## This has 8 parameters so its annoying to tune.

printf "FAS:\n" >> out5
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_order 2 -sub_fas_coarse_snes_max_it 1 -sub_fas_coarse_ksp_rtol 1e-1 -sub_fas_coarse_pc_type eisenstat >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_order 2 -sub_fas_coarse_snes_max_it 1 -sub_fas_coarse_ksp_rtol 1e-1 -sub_fas_coarse_pc_type eisenstat >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_order 2 -sub_fas_coarse_snes_max_it 1 -sub_fas_coarse_ksp_rtol 1e-1 -sub_fas_coarse_pc_type eisenstat >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_order 2 -sub_fas_coarse_snes_max_it 1 -sub_fas_coarse_ksp_rtol 1e-1 -sub_fas_coarse_pc_type eisenstat >> out5
##

printf "Anderson with LS NPC:\n" >> out5
lsnpc=' -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_order 2 -sub_npc_snes_max_it 1'
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc >> out5

printf "NGMRES with LS NPC:\n" >> out5
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc >> out5
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -sub_snes_max_it 1 -snes_converged_reason -snes_monitor -snes_view
## NGMRES doensn't converge without a preconditioner
## This is good

# printf "\nNGMRES with TR NPC:\n" >> out5
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_max_it 1 -sub_npc_snes_type newtontr -sub_npc_pc_type eisenstat -sub_npc_snes_max_it 1 format
## This is bad. Trustregion is not a good preconditioner.

printf "QN with LS NPC:\n" >> out5
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type qn -sub_snes_linesearch_order 2 -sub_snes_linesearch_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type qn -sub_snes_linesearch_order 2 -sub_snes_linesearch_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type qn -sub_snes_linesearch_order 2 -sub_snes_linesearch_max_it 1 $lsnpc >> out5
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_type qn -sub_snes_linesearch_order 2 -sub_snes_linesearch_max_it 1 $lsnpc >> out5
## It's not fast, but it has fewer DDM iterations.