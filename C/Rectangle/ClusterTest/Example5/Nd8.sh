#!/bin/bash -l
#SBATCH -J Nd8
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 16:0:0
#SBATCH --nodes 8
#SBATCH --ntasks 8
#SBATCH --ntasks-per-node 1
#SBATCH -o Nd8.slurm
#SBATCH -e Nd8.slurm
# SBATCH --exclusive

module load gnu8 mpich petsc/3.12.0
Nd=8
gtol=1e-6
ntol=1e-5
ktol=1e-3
op=0.4
h=0.01

file=Nd$Nd
rm -f $file.out

for limit in {1.5,2.0}
do
   N=$(echo "scale = 0; 2*$limit/$h" | bc)
   printf "now running [mpiexec -np $Nd ../../maddm -N $N -problem ex5 -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol $gtol]\n" >> $file.out
   mpiexec -np $Nd ../../maddm -htn -N $N -problem ex5 -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol $gtol >> $file.out
done