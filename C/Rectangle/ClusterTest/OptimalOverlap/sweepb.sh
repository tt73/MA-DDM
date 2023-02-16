#!/bin/bash -l
#SBATCH -J sweep-b
#SBATCH -o slurmb
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 24:0:0
#SBATCH --mem=0
#SBATCH --ntasks 4
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f sweepb.out

# N=200
np=4
prob=ex1
limit=2.0
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)

for ktol in {1e-2,1e-3,1e-4}
do
   for op in {0.025,0.05,0.075,0.10}
   do
      printf "Running: mpiexec -np $np ../../maddm -sin -N $N -op $op -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -sub_ksp_rtol $ktol -snes_converged_reason -snes_max_it 1000 \n" >> sweepb.out
      mpiexec -np $np ../../maddm -sin -N $N -op $op -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -sub_ksp_rtol $ktol -snes_converged_reason  -snes_max_it 1000 >> sweepb.out
   done
done