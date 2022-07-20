#!/bin/bash -l
## Just trying out different global solvers.

## NAMING
#SBATCH -J j6
#SBATCH -p public
#SBATCH -o slurmout6
#SBATCH -e slurmout6

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

rm -f out6

N=200
np=4

printf "NASM:\n" >> out6
mpiexec -np $np ../../maddm -N $N -problem ex1 -sin >> out6
mpiexec -np $np ../../maddm -N $N -problem ex2 -sin >> out6
mpiexec -np $np ../../maddm -N $N -problem ex3 -sin >> out6
mpiexec -np $np ../../maddm -N $N -problem ex4 -sin >> out6

printf "Newton-Krylov Schwarz:\n" >> out6
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex4 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2 >> out6

printf "ASPIN:\n" >> out6
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type aspin -npc_sub_pc_type ilu >> out6
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type aspin -npc_sub_pc_type ilu >> out6
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type aspin -npc_sub_pc_type ilu >> out6
mpiexec -np $np ../../maddm -N $N -problem ex4 -snes_type aspin -npc_sub_pc_type ilu >> out6

printf "FAS:\n" >> out6
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type fas -fas_coarse_snes_linesearch_order 2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type fas -fas_coarse_snes_linesearch_order 2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type fas -fas_coarse_snes_linesearch_order 2 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex4 -snes_type fas -fas_coarse_snes_linesearch_order 2 >> out6

printf "NGMRES with Newton NPC:\n" >> out6
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_order 2 -npc_snes_max_it 1 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_order 2 -npc_snes_max_it 1 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_order 2 -npc_snes_max_it 1 >> out6
mpiexec -np $np ../../maddm -N $N -problem ex4 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_order 2 -npc_snes_max_it 1 >> out6