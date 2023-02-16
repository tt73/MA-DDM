#!/bin/bash -l
#SBATCH -J sweep-a
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

rm -f sweepa.out

# N=200
np=4
prob=ex1
limit=1.5
h=0.01
N=$(echo "scale = 0; 2*$limit/$h" | bc)

for ntol in {1e-3,1e-4}
do
   for op in {0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40}
   do
      printf "Running: mpiexec -np $np ../../maddm -sin -N $N -op $op -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -sub_snes_rtol $ntol -sub_ksp_rtol 1e-1 -snes_converged_reason -snes_max_it 1000\n" >> sweepa.out
      mpiexec -np $np ../../maddm -htn -N $N -op $op -problem $prob -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -sub_snes_rtol $ntol -sub_ksp_rtol 1e-1 -snes_converged_reason -snes_max_it 1000 >> sweepa.out
   done
done