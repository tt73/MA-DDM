#!/bin/bash -l
#SBATCH -J NKS
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -o slurmout_nasm
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive

## Try to do NASM with local NKS
Nd=4
N=300
module load gnu8 mpich petsc/3.12.0


## Serial version
../../maddm -N $N -problem ex1

## NKS
mpiexec -np $Nd ../../maddm -N $N -nks -problem ex1 -snes_view


## NASM with usual NK local
mpiexec -np $Nd ../../maddm -N $N -HTN -problem ex1 -snes_view

## NASM with NKS locally
mpiexec -np $Nd ../../maddm -N $N -HTN -problem ex1 -snes_view -sub_pc_type asm -sub_pc_asm_blocks 2


