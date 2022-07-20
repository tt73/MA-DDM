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

N=200
np=4

printf "NASM:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sin
mpiexec -np $np ../../maddm -N $N -problem ex2 -sin
mpiexec -np $np ../../maddm -N $N -problem ex3 -sin
mpiexec -np $np ../../maddm -N $N -problem ex4 -sin

printf "Newton:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2

printf "ASPIN:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type aspin -npc_sub_pc_type ilu
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type aspin -npc_sub_pc_type ilu
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type aspin -npc_sub_pc_type ilu


printf "FAS:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type fas -fas_coarse_snes_linesearch_type l2

# printf "\nNGMRES:\n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres -snes_ngmres_select_type difference




printf "\nNGMRES with Newton LS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_type basic -npc_snes_max_it 1


printf "\nNGMRES with FAS NPC:\n"
fasnpc='-npc_snes_type fas -npc_fas_coarse_snes_linesearch_type l2 -npc_snes_max_it 1 -npc_fas_levels_snes_type newtonls -npc_fas_levels_snes_max_it 6'
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres $fasnpc


printf "\nFAS * NGMRES composite:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type composite -snes_composite_type multiplicative -snes_composite_sneses fas,ngmres

printf "\nFAS + Newton composite:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type composite -snes_composite_type additiveoptimal -snes_composite_sneses newtontr,newtonls