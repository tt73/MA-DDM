#!/bin/bash -l
#SBATCH -J Nd1
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd1.slurm
#SBATCH -e Nd1.slurm

Nd=1
module load gnu8 mpich petsc/3.12.0
rm -f test1.out

for limit in {0.5,1.0,1.5,2.0}
do
   printf "[-$limit, $limit]^2\n"
   for h in {0.05,0.01}
   do
      N=$(echo "scale = 0; 2*$limit/$h" | bc)
      mpiexec -np $Nd ../../maddm -N $N -problem ex1 -snes_rtol 1e-8 -ksp_rtol 1e-4 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> test1.out
   done
done

## I want to run the serial algorithm without any domain decomposition.
## The error should grow for fixed tol 1e-8 for newton and 1e-5 for gmres as the limits increase.