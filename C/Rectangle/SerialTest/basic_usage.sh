#!/bin/bash
# The purpose of this script is to show how the MA-DDM program can be run.
#
# Run this code with `bash basic_usage.sh`
# You can pass an integer as an arg to change N.


# N = Number of interior points
N=20

printf "Running code serially with $N^2 points on a square.\nNewton iterations: \n"
../maddm -t1_N ${1:-$N} -snes_monitor -snes_type newtonls
## Change the outer nonlinear solver with -snes_type.
## Turn on option to print residues with -snes_monitor.



printf "\n\nRunning code on rectangular domain [-.5, .5] x [-2, 2]:\n"
../maddm -t1_N ${1:-$N} -xmin -0.5 -xmax 0.5 -ymin -2.0 -ymax 2.0 -print
## The limits of the rectangular domain can be changed with



# Change problem
printf "\n\nTest3: Change the problem :\n"
printf "example 10: "; ../maddm -t1_N 10
printf "example 11: "; ../maddm -t1_N 10 -t1_problem ex11
printf "example 12: "; ../maddm -t1_N 10 -t1_problem ex12
## You can switch between 3 different convex problems.


## Change initial guess
printf "\n\nTest4: Try out initial different guesses for the Newton iteration :\n"
printf "zeros:   "; ../maddm -t1_N 40 -snes_converged_reason -snes_type newtontr -t1_init_type zeros | grep 'Nonlinear'
printf "random:  "; ../maddm -t1_N 40 -snes_converged_reason -snes_type newtontr -t1_init_type zeros | grep 'Nonlinear'
printf "corner:  "; ../maddm -t1_N 40 -snes_converged_reason -snes_type newtontr -t1_init_type corner | grep 'Nonlinear'
printf "pyramid: "; ../maddm -t1_N 40 -snes_converged_reason -snes_type newtontr -t1_init_type pyramid | grep 'Nonlinear'
## Zeros is the simplest, u0 = 0.
## Random sets the initial guess to be uniform random.
## Corner sets the u0 = M, where M is the largest value of g on the 4 corners.
## Pyramid sets the initial guess to be an inverted pyramid.


# Convergence test
printf "\n\nTest 5: Convergence wrt N\n"
for i in {10..50..10}
do
   ../maddm -t1_N $i -snes_type fas
done
## the error descreases as expected
## the epsilon and width dynamically with N


printf "\n\nTest 6: Various SNES types\n"
printf "List of availables SNES types:\n"
../maddm -help | grep '-snes_type '
snes_types="newtonls newtontr nrichardson ksponly ksptransposeonly vinewtonssls ngmres qn ncg fas nasm anderson aspin"
for s in $snes_types
do
   printf "Trying out ${s}: "
   ../maddm -t1_N ${1:-$N} -snes_type $s
done
## You can see which methods worked based on the error.
## Converge: newtonls newtontr vinewtonssls fas nasm aspin
## Diverge: nrichardson ksponly ksptransposeonly ngmres qn ncg anderson
##
## ksponly diverges because that is only meant for linear problems.
## nrichardson, ngmres, qn, ncg, and anderson all converge for the bratue example. Not sure whey it doesn't work for our code.
## nasm and aspin are parallel codes so they are effectively the same as newtonls in the serial case.


printf "\n\nTest 7: Various SNES types with Nonlinear preconditioning\n"
snes_types="newtonls newtontr nrichardson ksponly ksptransposeonly vinewtonssls ngmres qn ncg nasm anderson"
for s in $snes_types
do
   printf "Trying out ${s}: "
   ../maddm -t1_N ${1:-$N} -snes_type $s -npc_snes_type fas
done
## The snes types that did not converge earlier: nrichardson, ngmres, qn, ncg, and anderson all converge.
## The



printf "\n\nTest 8: Various KSP types\n"
printf "List of availables SNES types:\n"
../maddm -help | grep 'ksp_type'
ksp_types="cg groppcg pipecg pipecgrr pipelcg pipeprcg pipecg2 cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls"
for k in $ksp_types
do
   printf "Trying out $k: "
   ../maddm -t1_N ${1:-$N} -snes_type newtonls -ksp_type $k | grep '*Error'
done
## We use newton linesearch for the nonlinear solver, and vary the linear solver option.
## Converge: cgne richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr preonly bicg fgmres pipefgmres lgmres lcd gcr pipegcr pgmres dgmres cgls


printf "\n\nTest 9: Various PC types\n"
printf "List of availables PC types:\n"
../maddm -help | grep 'pc_type'
pc_types="none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat cp redistribute svd gamg kaczmarz hmg lmvm deflation"
for p in $pc_types
do
   printf "Trying out $p: "
   ../maddm -t1_N ${1:-$N} -snes_type newtonls -pc_type $p | grep '*Error'
done
## We use newton linesearch for the nonlinear solver, gmres for the Krylov method, and vary the preconditioner.
## Converge: none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat redistribute svd gamg kaczmarz hmg deflation