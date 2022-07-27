#!/bin/bash -l
## Single iteration newton
## We keep newton iteration at 1, and vary the GMRES tolerance.

## NAMING
#SBATCH -J j2
#SBATCH -o slurmout2
#SBATCH -e slurmout2

## partition/queue
#SBATCH -p dms-cpu

## EMAIL NOTIFICATION
#SBATCH --mail-user tt73@njit.edu
#SBATCH --mail-type=END
#SBATCH -A tt73

## RUNTIME HOURS:MIN:SEC and MEMORY
#SBATCH -t 8:0:0
#SBATCH --mem=0G

## Task allocation
#SBATCH --ntasks 4
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f out2

N=200
np=4

for ksptol in {1.e-4,1.e-3,1.e-2,1.e-1}
do
   printf "ksptol = $ksptol\n" >> out2
   mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol >> out2
   mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol >> out2
   mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol >> out2
   mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol >> out2
done
