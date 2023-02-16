#!/bin/bash -l
#SBATCH -J NKS
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive

## This is a tutorial for the NKS with additive schwarz
Nd=4
N=200
module load gnu8 mpich petsc/3.12.0
rm -f NKS.out


## usual block jacobi, do not specify a PC
mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -ksp_type pipefgmres -snes_max_it 1000 -snes_view -snes_monitor

## additive schwarz PC
# mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -ksp_type pipefgmres -pc_type asm -snes_max_it 1000 -snes_view -snes_monitor >> NKS.out


## asm options
# ../../maddm -nks -snes_type newtonls -ksp_type pipefgmres -pc_type asm -help | grep 'pc_asm'>> NKS.out

# printf "\n\nTesting overlap for restrictive additive schwarz...\n" >> NKS.out
# for op in {0,5,10,15,20,25,30}
# do
#    printf "overlap = $op\n" >> NKS.out
#    mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -ksp_type pipefgmres -pc_type asm -snes_max_it 1000 -pc_asm_overlap $op >> NKS.out
# done

# printf "\n\nTesting overlap for basic additive schwarz...\n" >> NKS.out
# for op in {0,5,10,15,20,25,30}
# do
#    printf "overlap = $op\n" >> NKS.out
#    mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -ksp_type pipefgmres -pc_type asm -snes_max_it 1000 -pc_asm_overlap $op -pc_asm_type basic >> NKS.out
# done

# printf "\n\nTesting overlap for interpolating additive schwarz...\n" >> NKS.out
# for op in {0,5,10,15,20,25,30}
# do
#    printf "overlap = $op\n" >> NKS.out
#    mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_type newtonls -ksp_type pipefgmres -pc_type asm -snes_max_it 1000 -pc_asm_overlap $op -pc_asm_type interpolate >> NKS.out
# done



