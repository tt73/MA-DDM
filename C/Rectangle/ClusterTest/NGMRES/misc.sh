#!/bin/bash -l
## Test on 4 subdomains
N=300

#SBATCH -J misc
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 8:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
module load gnu8 mpich petsc/3.12.0

rm -f misc

## NGMRES -L SIN
mpirun ../../maddm -N $N -problem ex1 -ngmres -snes_npc_side left >> misc
mpirun ../../maddm -N $N -problem ex2 -ngmres -snes_npc_side left >> misc
mpirun ../../maddm -N $N -problem ex3 -ngmres -snes_npc_side left >> misc
mpirun ../../maddm -N $N -problem ex4 -ngmres -snes_npc_side left >> misc

## NGMRES -R SIN
mpirun ../../maddm -N $N -problem ex1 -ngmres -snes_npc_side right >> misc
mpirun ../../maddm -N $N -problem ex2 -ngmres -snes_npc_side right >> misc
mpirun ../../maddm -N $N -problem ex3 -ngmres -snes_npc_side right >> misc
mpirun ../../maddm -N $N -problem ex4 -ngmres -snes_npc_side right >> misc

## FAS + SIN
mpirun ../../maddm -N $N -problem ex1 -coarse >> misc
mpirun ../../maddm -N $N -problem ex2 -coarse >> misc
mpirun ../../maddm -N $N -problem ex3 -coarse >> misc
mpirun ../../maddm -N $N -problem ex4 -coarse >> misc