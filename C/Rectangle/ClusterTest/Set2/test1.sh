#!/bin/bash -l
## This is the initial test to make the local solver faster.
## Want to find out how far we can relax the newton and kyrlov tolerances before the method suffers.

## NAMING
#SBATCH -J j1
#SBATCH -o slurmout1
#SBATCH -e slurmout1

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

module load gnu8 mpich petsc

rm -f out1

N=200
np=4

for snestol in {1.e-4,1.e-3,1.e-2,1.e-1}
do
   for ksptol in {1.e-4,1.e-3,1.e-2,1.e-1}
   do
      printf "snestol = $snestol, ksptol = $ksptol\n" >> out1
      mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol >> out1
      mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol >> out1
      mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol >> out1
      mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol >> out1
   done
done