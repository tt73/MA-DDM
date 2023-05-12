#!/bin/bash -l
#SBATCH -J ex1-iters
#SBATCH -p lowpriority
#SBATCH --constraint=avx512
#SBATCH --mail-user tt73@njit.edu
#SBATCH --mail-type=END
#SBATCH -A tt73
#SBATCH -t 24:0:0
#SBATCH --mem=0G
#SBATCH -N 1

## These are the tests for Example 1 on [-L,L]^2
## The main parameters are:
##  * Domain limits: L = 0.5, 1.0, 1.5, 2.0
##  * Overlap percentage: p = 0.1, 0.2, 0.3, 0.4
##  * Mesh resolution: h = 0.05, 0.01
##  * Decomposition in n by m format: 2x1, 2x2, 3x2, 4x2, 3x3, 4x4, 5x5, 6x6

## The initial guess for all runs is the coarse solve on a grid with resolution 4h
## If we truly wanted to solve resolution h, then we must choose N = 2*L/h-1,
## but for the compatitibility for the initilization we have to choose 2*L/h+1 (bigger by 2).

module load gnu8 mpich petsc/3.12.0

rm -f out

for Nd in {2,4,6,8,9,16,25,36}
do
   for op in {0.1,0.2,0.3,0.4}
   do
      for limit in {0.5,1.0,1.5,2.0}
      do
         for h in {0.05,0.01}
         do
            N=$(echo "scale = 0; 2*$limit/$h+1" | bc)
            mpiexec -np $Nd ../maddm -N $N -problem ex1 -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -op $op -init_type coarse -coarseness 4 -snes_converged_reason >> out
         done
      done
   done
done