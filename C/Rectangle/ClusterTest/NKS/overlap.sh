#!/bin/bash -l
#SBATCH -J overlap
#SBATCH -o slurmout1
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 8:0:0
#SBATCH --mem=0
#SBATCH --ntasks 4
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f overlap.out

N=200
np=4
for op in {0.0,0.05,0.10,0.15,0.20,0.25,0.30}
do
   printf "overlap = $op\n" >> overlap.out
   mpirun ../../maddm -N $N -op $op -problem ex1 -nks >> overlap.out
   mpirun ../../maddm -N $N -op $op -problem ex2 -nks >> overlap.out
   mpirun ../../maddm -N $N -op $op -problem ex3 -nks >> overlap.out
   mpirun ../../maddm -N $N -op $op -problem ex4 -nks >> overlap.out
done