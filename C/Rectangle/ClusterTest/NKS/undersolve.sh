#!/bin/bash -l
#SBATCH -J undersolve
#SBATCH -p dms-cpu
#SBATCH -o undersolve.slurm
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive

## Relax the tolerance on the Krylov
Nd=4
N=300
module load gnu8 mpich petsc/3.12.0
rm -f undersolve.out


for ksptol in {1.e-5,1.e-4,1.e-3,1.e-2,1.e-1}
do
   printf "ksptol = $ksptol\n" >> undersolve.out
   mpiexec -np $Nd ../../maddm -N $N -problem ex1 -nks -ksp_rtol $ksptol -op 0 >> undersolve.out
   mpiexec -np $Nd ../../maddm -N $N -problem ex2 -nks -ksp_rtol $ksptol -op 0 >> undersolve.out
   mpiexec -np $Nd ../../maddm -N $N -problem ex3 -nks -ksp_rtol $ksptol -op 0 >> undersolve.out
   mpiexec -np $Nd ../../maddm -N $N -problem ex4 -nks -ksp_rtol $ksptol -op 0 >> undersolve.out
done
