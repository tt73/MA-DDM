N=200
np=4

## Linesearch types ============================
## shell basic l2 bt nleqerr cp ncglinear

printf "\nlinesearch bt - 2nd order:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## slightly faster than 3rd order method

printf "\nlinesearch l2:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## I tried changing the linesearch method to a lower order method.
## Number of iterations doesn not change.
## It's actually a bit faster than bt used in the base settings in terms of wtime.

printf "\nlinesearch basic:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Takes full step every time
## It dosen't converge unless you have the damping < 1.

printf "\nInexact method: \n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_ksp_type preonly -sub_pc_type ilu -sub_snes_linesearch_type bt -sub_snes_max_it 1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_ksp_type preonly -sub_pc_type ilu -sub_snes_linesearch_type bt -sub_snes_max_it 1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_ksp_type preonly -sub_pc_type ilu -sub_snes_linesearch_type bt -sub_snes_max_it 1
## Works suprisingly well for problem 3