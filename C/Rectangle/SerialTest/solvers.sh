#!/bin/bash
# The purpose of this script is to show how the MA-DDM program can be run.
#
# Run this code with `bash basic_usage.sh`
# You can pass an integer as an arg to change N.


# N = Number of interior points
N=20

printf "\n\nVarious SNES types\n"
printf "List of availables SNES types:\n"
../maddm -help | grep '-snes_type '
snes_types="newtonls newtontr nrichardson vinewtonssls ngmres qn ncg fas nasm anderson aspin"
for s in $snes_types
do
   printf "Trying out ${s}: "
   ../maddm -N ${1:-$N} -snes_type $s
done
## You can see which methods worked based on the error.
## Converge: newtonls newtontr vinewtonssls fas nasm aspin
## Diverge: nrichardson ngmres qn ncg anderson
##
## nrichardson, ngmres, qn, ncg, and anderson all converge for the bratue example. Not sure whey it doesn't work for our code.
## nasm and aspin are parallel codes so they are effectively the same as newtonls in the serial case.


printf "Various SNES types with Nonlinear preconditioning\n"
snes_types="newtonls newtontr nrichardson ksponly ksptransposeonly vinewtonssls ngmres qn ncg nasm anderson"
for s in $snes_types
do
   printf "Trying out ${s}: "
   ../maddm -N ${1:-$N} -snes_type $s -npc_snes_type fas
done
## The snes types that did not converge earlier: nrichardson, ngmres, qn, ncg, and anderson all converge.


printf "\n\nVarious KSP types\n"
printf "List of availables SNES types:\n"
../maddm -help | grep 'ksp_type'
ksp_types="cg groppcg pipecg pipecgrr pipelcg pipeprcg pipecg2 cgne nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls"
for k in $ksp_types
do
   printf "Trying out $k: "
   ../maddm -N ${1:-$N} -snes_type newtonls -ksp_type $k | grep '*Error'
done
## We use newton linesearch for the nonlinear solver, and vary the linear solver option.
## Converge: cgne richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr preonly bicg fgmres pipefgmres lgmres lcd gcr pipegcr pgmres dgmres cgls


printf "\n\nVarious PC types\n"
printf "List of availables PC types:\n"
../maddm -help | grep 'pc_type'
pc_types="none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat cp redistribute svd gamg kaczmarz hmg lmvm deflation"
for p in $pc_types
do
   printf "Trying out $p: "
   ../maddm -N ${1:-$N} -snes_type newtonls -pc_type $p | grep '*Error'
done
## We use newton linesearch for the nonlinear solver, gmres for the Krylov method, and vary the preconditioner.
## Converge: none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat redistribute svd gamg kaczmarz hmg deflation