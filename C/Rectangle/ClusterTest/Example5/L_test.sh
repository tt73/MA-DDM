#!/bin/bash -l
#SBATCH -J L
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 24:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o slurm_L
#SBATCH -e slurm_L

## Test on 4 subdomains
Nd=4
module load gnu8 mpich petsc/3.12.0
rm -f out_L

for limit in {0.5,1.0,1.5,2.0}
do
   for h in {0.05,0.01}
   do
      N=$(echo "scale = 0; 2*$limit/$h" | bc)
      mpiexec -np $Nd ../../maddm -N $N -problem ex5 -htn -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-8 >> out_L
      mpiexec -np $Nd ../../maddm -N $N -problem ex5 -sin -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-8 >> out_L
      mpiexec -np $Nd ../../maddm -N $N -problem ex5 -nks -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-8 >> out_L
   done

done