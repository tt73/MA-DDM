#!/bin/bash -l
#SBATCH -J Nd9
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 16:0:0
#SBATCH --nodes 9
#SBATCH --ntasks 9
#SBATCH --ntasks-per-node 1
#SBATCH -o Nd9.slurm
#SBATCH -e Nd9.slurm
# SBATCH --exclusive

module load gnu8 mpich petsc/3.12.0
Nd=9
gtol=1e-8
ntol=1e-3
ktol=1e-2
# op=0.4
h=0.01

file=Nd$Nd
rm -f $file.out

for op in {0.1,0.2,0.4}
do
   for limit in {1.5,2.0}
   do
      N=$(echo "scale = 0; 2*$limit/$h" | bc)
      printf "now running [mpiexec -np $Nd ../../maddm -htn -N $N -problem ex5 -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol $gtol]\n" >> $file.out
      mpiexec -np $Nd ../../maddm -htn -N $N -problem ex5 -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol $gtol >> $file.out
   done
done