#!/bin/bash -l
#SBATCH -J serial
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -o slurmout_serial
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive

module load gnu8 mpich petsc/3.12.0

## Serial version
for N in {100,150,200,250,300}
do
   printf "N  = $N\n"
   ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -pc_type eisenstat
   ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -pc_type asm -pc_asm_blocks 2
   ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -pc_type asm -pc_asm_blocks 4
   printf "\n"
done

