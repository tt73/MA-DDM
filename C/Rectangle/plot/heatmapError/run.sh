#!/bin/bash -l
#SBATCH -J generatesol
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 16:0:0
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --exclusive


module load gnu8 mpich petsc/3.12.0

gtol=1e-8
ntol=1e-3
ktol=1e-2
op=0.2
h=0.01
limit=2.0
N=$(echo "scale = 0; 2*$limit/$h" | bc)

for Nd in {2,4,6,8,9}
do
   mpiexec -np $Nd ../../maddm -htn -N $N -problem ex5 -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason -snes_rtol $gtol -sol
   mv load_u.m load_u$Nd.m
   mv load_exact.m load_exact$Nd.m
done
