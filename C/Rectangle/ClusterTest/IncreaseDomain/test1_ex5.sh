#!/bin/bash -l
#SBATCH -J test1ex5
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
rm -f test1ex5.out

for limit in {0.5,1.0,1.5,2.0}
do
   for h in {0.05,0.01}
   do
      N=$(echo "scale = 0; 2*$limit/$h" | bc)
      # echo "../../maddm -N $N -problem ex5 -snes_rtol 1e-8 -ksp_rtol 1e-5 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason "
      ../../maddm -N $N -problem ex5 -snes_rtol 1e-8 -ksp_rtol 1e-4 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> test1ex5.out
   done
done
