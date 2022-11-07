#!/bin/bash -l
#SBATCH -J Nd1
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 12:0:0
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd1.slurm
#SBATCH -e Nd1.slurm

## Test on 1 subdomain (serial)
Nd=1
module load gnu8 mpich petsc/3.12.0
rm -f Nd1.out

for limit in {0.5,1.0,1.5,2.0}
do
   printf "[-$limit, $limit]^2\n"
   for tol in {1e-1,1e-2,1e-3,1e-4,1e-5,1e-6}
   do
      printf "tol = $tol\n"
      for h in {0.05,0.01}
      do
         N=$(echo "scale = 0; 2*$limit/$h" | bc)  # it works
         printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex2c -sub_snes_rtol $tol -sub_ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit\n"
         ../../maddm -N $N -problem ex2c -snes_rtol $tol -ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> Nd1.out
      done
   done
done

printf "all done."