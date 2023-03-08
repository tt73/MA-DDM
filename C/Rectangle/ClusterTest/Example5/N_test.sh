#!/bin/bash -l
#SBATCH -J N
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 10:0:0
#SBATCH --nodes 1
#SBATCH -o slurm_N
#SBATCH -e slurm_N

rm -f out_N

for N in 100 150 200 250 300 350 400 450 500
do
   printf "N = $N\n" >> out_N
   ../../maddm -N $N -problem ex5 >> out_N
done