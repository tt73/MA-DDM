#!/bin/bash -l
#SBATCH -J Nd1
#SBATCH -p dms-cpu
#SBATCH --mail-user tt73@njit.edu
#SBATCH -A tt73
#SBATCH -t 8:0:0
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
   for tol in {1e-1,1e-4,1e-6}
   do
      printf "tol = $tol\n"
      for N in {100,200,300,400}
      do
         printf "Now running mpiexec -np $Nd ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit\n"
         ../../maddm -N $N -problem ex1 -sub_snes_rtol $tol -sub_ksp_rtol $tol -xmin -$limit -xmax $limit -ymin -$limit -ymax $limit >> Nd1.out
      done
   done
done

printf "all done."


# ## NGMRES -L SIN
# mpirun ../../maddm -N $N -problem ex1 -ngmres -snes_npc_side left >> misc
# mpirun ../../maddm -N $N -problem ex2 -ngmres -snes_npc_side left >> misc
# mpirun ../../maddm -N $N -problem ex3 -ngmres -snes_npc_side left >> misc
# mpirun ../../maddm -N $N -problem ex4 -ngmres -snes_npc_side left >> misc

# ## NGMRES -R SIN
# mpirun ../../maddm -N $N -problem ex1 -ngmres -snes_npc_side right >> misc
# mpirun ../../maddm -N $N -problem ex2 -ngmres -snes_npc_side right >> misc
# mpirun ../../maddm -N $N -problem ex3 -ngmres -snes_npc_side right >> misc
# mpirun ../../maddm -N $N -problem ex4 -ngmres -snes_npc_side right >> misc

# ## FAS + SIN
# mpirun ../../maddm -N $N -problem ex1 -coarse >> misc
# mpirun ../../maddm -N $N -problem ex2 -coarse >> misc
# mpirun ../../maddm -N $N -problem ex3 -coarse >> misc
# mpirun ../../maddm -N $N -problem ex4 -coarse >> misc