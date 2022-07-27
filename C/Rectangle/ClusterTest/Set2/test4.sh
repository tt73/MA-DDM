#!/bin/bash -l
## I'm trying out various linear solvers and preconditioners.
## It's all with HTN.

## NAMING
#SBATCH -J j4
#SBATCH -o slurmout4
#SBATCH -e slurmout4

## partition/queue
#SBATCH -p dms-cpu

## EMAIL NOTIFICATION
#SBATCH --mail-user tt73@njit.edu
#SBATCH --mail-type=END
#SBATCH -A tt73

## RUNTIME HOURS:MIN:SEC and MEMORY
#SBATCH -t 8:0:0
#SBATCH --mem=0G

## Task allocation
#SBATCH --ntasks 4
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 1

module load gnu8 mpich petsc/3.12.0

rm -f out4

N=200
np=4

printf "base of comparison - DGMRES + eisenstat:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4  -sin >> out4
## Tests 1, 2, and 3 were run with this type of subsolve.
## It's performs well.

# PC TYPES ============================
# none jacobi bjacobi sor lu mg eisenstat ilu asm gasm gamg
# I removed ones that don't work with GMRES from the list, like icc.

printf "GMRES + multigrid:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type mg  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type mg  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type mg  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type gmres -sub_pc_type mg  -sin >> out4
## Not that good.

printf "GMRES + jacobi:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type jacobi  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type jacobi  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type jacobi  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type gmres -sub_pc_type jacobi  -sin >> out4
## This one works buts its a bit slower.

printf "GMRES + SOR:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type sor  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type sor  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type sor  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type gmres -sub_pc_type sor  -sin >> out4
## This one is decent.


# KSP TYPES ============================
# richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls
#

printf "PFGMRES:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type pipefgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type pipefgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type pipefgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type pipefgmres  -sin >> out4
## pipelined flexible GMRES
## its close to GMRES speed

printf "Chebyshev:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type chebyshev  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type chebyshev  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type chebyshev  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type chebyshev  -sin >> out4
## Doesn't seem to work for problem 3
##

printf "GCR:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gcr  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gcr  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gcr  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type gcr  -sin >> out4
##

printf "LGMRES:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type lgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type lgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type lgmres  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type lgmres  -sin >> out4

printf "BCGSL:\n" >> out4
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type bcgsl  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type bcgsl  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type bcgsl  -sin >> out4
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_ksp_type bcgsl  -sin >> out4