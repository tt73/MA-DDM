#!/bin/bash -l
#SBATCH -J sweep-c
#SBATCH -o slurma
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 24:0:0
#SBATCH --mem=0
#SBATCH --ntasks 9
#SBATCH --nodes 9
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f sweepc.out

# N=200
np=9
prob=ex3
ntol=1e-1
ktol=1e-1

for N in {100,200,300,400}
do
   for op in {0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40}
   do
      printf "Running: mpiexec -np $np ../../maddm -htn -N $N -op $op -problem $prob -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -snes_converged_reason -snes_max_it 1000\n" >> sweepc.out
      mpiexec -np $np ../../maddm -htn -N $N -op $op -problem $prob -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -snes_converged_reason -snes_max_it 1000 >> sweepc.out
   done
done