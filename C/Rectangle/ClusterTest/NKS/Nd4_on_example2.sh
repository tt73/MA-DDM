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

## Test on 4 subdomains for example 1
## We want to see how NKS performs as the problem size gets bigger
Nd=4
module load gnu8 mpich petsc/3.12.0
rm -f Nd4.out

for limit in {0.5,1.0,1.5,2.0}
do
   printf "[-$limit, $limit]^2\n"
   for tol in 1e-6
   do
      printf "tol = $tol\n"
      for h in {0.05,0.01}
      do
         N=$(echo "scale = 0; 2*$limit/$h" | bc)
         printf "Now running mpiexec -np $Nd ../../maddm  -nks -N $N -problem ex2c -snes_type newtonls -snes_linesearch_order 2 -ksp_type pipefgmres -pc_type bjacobi -ksp_rtol $tol -snes_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit\n"
         mpiexec -np $Nd ../../maddm -N $N -nks  -problem ex2c -snes_type newtonls -snes_linesearch_order 2 -ksp_type pipefgmres -pc_type bjacobi -ksp_rtol $tol -snes_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_max_it 2000 -snes_converged_reason >> Nd4.out
      done
   done
done