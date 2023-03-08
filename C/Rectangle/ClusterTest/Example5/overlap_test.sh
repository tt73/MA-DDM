#!/bin/bash -l
#SBATCH -J overlap
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o slurm_overap
#SBATCH -e slurm_overap

## Test on 4 subdomains
Nd=4
module load gnu8 mpich petsc/3.12.0
rm -f out_overlap

N=200
np=4
for op in {0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45}
do
   printf "overlap = $op\n" >> out_overlap
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex5 -sin >> out_overlap
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex5 -htn >> out_overlap
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex5 -nks >> out_overlap
done