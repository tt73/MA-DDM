#!/bin/bash -l
#SBATCH -J sweep-NKS
#SBATCH -o slurma
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 24:0:0
#SBATCH --mem=0
#SBATCH --ntasks 4
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f sweepNKS.out

# N=200
np=4
prob=ex5
limit=2.0
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)

for ktol in {1e-1,1e-2,1e-3,1e-4,1e-5}
do
   printf "Running: mpiexec -np $np ../../maddm -nks -N $N -op 0 -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -sub_ksp_rtol $ktol -snes_converged_reason -snes_max_it 1000\n" >> sweepNKS.out
   mpiexec -np $np ../../maddm -nks -N $N -op 0 -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -ksp_rtol $ktol -snes_converged_reason -snes_max_it 1000 >> sweepNKS.out

done