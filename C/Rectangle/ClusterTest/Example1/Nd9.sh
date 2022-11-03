#!/bin/bash -l
#SBATCH -J Nd9
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 16:0:0
#SBATCH --nodes 9
#SBATCH --ntasks 9
#SBATCH --ntasks-per-node 1
#SBATCH --exclusive
#SBATCH -o Nd9.slurm
#SBATCH -e Nd9.slurm

## Test on 4 subdomains
module load gnu8 mpich petsc/3.12.0
Nd=9
rm -f Nd9.out

for limit in {1.5,2.0}
do
   printf "[-$limit, $limit]^2\n"
   for tol in 1e-6
   do
      printf "tol = $tol\n"
      for h in {0.05,0.01}
      do
         N=$(echo "scale = 0; 2*$limit/$h" | bc)
         printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason\n"
         mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit -snes_converged_reason >> Nd9.out
      done
   done
done

printf "done"