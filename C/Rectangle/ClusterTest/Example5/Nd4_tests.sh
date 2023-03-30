#!/bin/bash -l
#SBATCH -J Nd4
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd4_tests.slurm
#SBATCH -e Nd4_tests.slurm

## Test on 4 subdomains
Nd=4
module load gnu8 mpich petsc/3.12.0
rm -f Nd4_tests.out

limit=2.0
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)
tol=1e-1

printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason\n"
mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> Nd4_test.out


printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-10\n"
mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-10 >> Nd4_test.out
