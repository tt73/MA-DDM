#!/bin/bash -l
## This is the initial test to make the local solver faster.
## This is mostly about testing the linesearch methods.
## 3a is for HTN and 3b is for SIN.

## NAMING
#SBATCH -J j3
#SBATCH -p public
#SBATCH -o slurmout3
#SBATCH -e slurmout3

## partition/queue
#SBATCH -p dms-cpu

## EMAIL NOTIFICATION
#SBATCH --mail-user tt73@njit.edu
#SBATCH --mail-type=END
#SBATCH -A tt73

## RUNTIME HOURS:MIN:SEC and MEMORY
#SBATCH -t 8:0:0
#SBATCH --mem=16G
#SBATCH -N 4

N=200
np=4

rm -f out3a out3b

## Linesearch types ============================
## shell basic l2 bt nleqerr cp ncglinear

printf "linesearch bt - 2nd order:\n" | tee --append out3a out3b
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a

mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_linesearch_type bt -sub_snes_linesearch_order 2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
## slightly faster than 3rd order method

printf "linesearch l2:\n" | tee --append out3a out3b
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_linesearch_type l2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a

mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex4 -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
## I tried changing the linesearch method to a lower order method.
## Number of iterations doesn not change.
## It's actually a bit faster than bt used in the base settings in terms of wtime.

printf "linesearch basic:\n" | tee --append out3a out3b
mpiexec -np $np ../../maddm -N $N -problem ex1  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex2  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex3  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a
mpiexec -np $np ../../maddm -N $N -problem ex4  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 >> out3a

mpiexec -np $np ../../maddm -N $N -problem ex1  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex2  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex3  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
mpiexec -np $np ../../maddm -N $N -problem ex4  -sub_snes_linesearch_type basic -sub_snes_linesearch_damping 0.5 -sub_snes_max_it 1 -sub_ksp_rtol 1e-1 >> out3b
## Takes full step every time.
## It dosen't converge unless you have the damping < 1.