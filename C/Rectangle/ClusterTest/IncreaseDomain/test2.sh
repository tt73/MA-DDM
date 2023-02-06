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
#SBATCH -o Nd4.slurm
#SBATCH -e Nd4.slurm

## Test on 4 subdomains
Nd=4
module load gnu8 mpich petsc/3.12.0
rm -f test2.out

for limit in {0.5,1.0,1.5,2.0}
do
   printf "[-$limit, $limit]^2\n"
   for h in {0.05,0.01}
   do
      N=$(echo "scale = 0; 2*$limit/$h" | bc)
      mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sin -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol 1e-8 >> test2.out
   done

done