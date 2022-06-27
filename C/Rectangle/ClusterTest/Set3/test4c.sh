N=200
np=4

printf "base of comparison - DGMRES + eisenstat:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Tests 1, 2, and 3 were run with this type of subsolve.
## It's performs well.

# PC TYPES ============================
# none jacobi bjacobi sor lu mg eisenstat ilu asm gasm gamg
# I removed ones that don't work with GMRES from the list, like icc.


printf "\nGMRES + multigrid:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type mg  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type mg  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type mg  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Not that good.

printf "\nGMRES + jacobi:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type jacobi  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type jacobi  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type jacobi  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## This one works buts its a bit slower.

printf "\nGMRES + SOR:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gmres -sub_pc_type sor  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gmres -sub_pc_type sor  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gmres -sub_pc_type sor  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## This one is decent.



# KSP TYPES ============================
# richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls
#

printf "\nPFGMRES:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type pipefgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type pipefgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type pipefgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## pipelined flexible GMRES
## its close to GMRES speed

printf "\nChebyshev:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type chebyshev  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type chebyshev  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type chebyshev  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Doesn't seem to work for problem 3
##

printf "\nGCR:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type gcr  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type gcr  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type gcr  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
##

printf "\nLGMRES:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type lgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type lgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type lgmres  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1

printf "\nBCGSL:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type bcgsl  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type bcgsl  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type bcgsl  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1