#!/bin/bash -l
#SBATCH -J Nd6
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 6
#SBATCH --ntasks 6
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd6_tests.slurm
#SBATCH -e Nd6_tests.slurm

## Test on 6 subdomains
Nd=6
module load gnu8 mpich petsc/3.12.0
rm -f Nd6_tests.out

limit=1.0
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)
tol=1e-1

printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_monitor\n"
mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_monitor >> Nd6_tests.out


# printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-10\n"
# mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol 1e-1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-10 >> Nd6_tests.out
