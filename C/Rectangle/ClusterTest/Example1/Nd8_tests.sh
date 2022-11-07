#!/bin/bash -l
#SBATCH -J Nd8
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 8
#SBATCH --ntasks 8
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd8_tests.slurm
#SBATCH -e Nd8_tests.slurm

Nd=8
module load gnu8 mpich petsc/3.12.0
rm -f Nd8_tests.out

limit=2.0
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)
tol=1e-2

# printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason\n"
# mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> Nd8_tests.out


printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-9\n"
mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-9 -snes_monitor >> Nd8_tests.out
