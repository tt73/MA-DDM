#!/bin/bash -l
#SBATCH -J j4
#SBATCH -o slurmout4
#SBATCH -p dms-cpu
#SBATCH -A tt73
#SBATCH -t 8:0:0
#SBATCH --mem=0G
#SBATCH --nodes 4
#SBATCH --ntasks 4
#SBATCH --ntasks-per-node 1
module load gnu8 mpich petsc/3.12.0

## the goal is to show that the iterating with Newton toleracnce down to tol h is not optimal
##
rm -f out
problem=ex1
N=201
L=1.0
op=0.2
h=$(echo "scale = 6; 2*$L/($N+1)" | bc)
## first, do a serial run
printf "serial with h = $h\n" >> out
../../../Rectangle/maddm -N $N -problem $problem -xmin -$L -xmax $L -ymin -$L -ymax $L -init_type coarse -coarseness 4 -snes_converged_reason -op $op  >> out
   printf "done\n\n\n" >> out


## parallel run
printf "parallel\n" >> out
for i in {1,2,4,8,16,32,64,128,256,512,1024,2048}
do
   tol=$(echo "scale = 6; $h * $i" | bc)
   printf "tol = $tol\n" >> out
   mpirun ../../../Rectangle/maddm -N $N -problem $problem -xmin -$L -xmax $L -ymin -$L -ymax $L -init_type coarse -coarseness 4 -snes_converged_reason -op $op -sub_snes_atol $tol -sub_snes_converged_reason >> out
    printf "done\n\n\n" >> out
done