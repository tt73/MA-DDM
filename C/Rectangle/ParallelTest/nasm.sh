#!/bin/bash

# N = Number of interior points
N=20

# We want to solve det(D^2u)=f using schwarz on the outside to decompose
# the domain and then use some non-linear solver on the local problems.
# The partitioning of the domain is automatic. Combining of the local
# solutions is also automatic. It is handled with snes type NASM.
printf "\n\nNon-linear additive schwarz options\n"
../test1 -help | grep nasm
#   -snes_nasm_type <now BASIC : formerly BASIC> Type of restriction/extension (choose one of) NONE RESTRICT INTERPOLATE BASIC ()
#   -snes_nasm_damping <1. : 1.>: The new solution is obtained as old solution plus dmp times (sum of the solutions on the subdomains) (SNESNASMSetDamping)
#   -snes_nasm_finaljacobian: <FALSE : FALSE> Compute the global jacobian of the final iterate (for ASPIN) ()
#   -snes_nasm_finaljacobian_type <now FINALOUTER : formerly FINALOUTER> The type of the final jacobian computed. (choose one of) FINALOUTER FINALINNER INITIAL ()
#   -snes_nasm_log: <TRUE : TRUE> Log times for subSNES solves and restriction ()


printf "\n\nTest 1: No DDM, just linesearch\n"
mpiexec -n 1 ../test1 -t1_N ${1:-$N} -snes_view
mpiexec -n 1 ../test1 -t1_N ${1:-$N} -log_view | grep 'problem\|error\|Nonlinear\|Time (sec):'
## The -snes_view shows exactly what's happening on the subdomains. They have the prefix
## Use the output for this as reference for the error

printf "\n\nTest 2: Two domains\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -da_overlap 0
## It converges. For N = 20, it takes 20 schwarz iterations to converge. We can set overlap to 0 (whatever that means).


printf "\n\nTest 3: Two domains again with simple single-step newton for each subdomain\n"
mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -sub_snes_linesearch_type basic -sub_snes_linesearch_max_it 1
printf "Default options   : "; mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_type nasm -log_view | grep 'Time (sec):'
printf "Single full newton: "; mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_type nasm -sub_snes_linesearch_type basic -sub_snes_linesearch_max_it 1 -log_view | grep 'Time (sec):'
## Instead of doing the cubic linesearch newton, we can opt to do a single full-step newton for the local problems. The speed difference is negligible.
## Perhaps I have to reconfigure PETSc with debugging mode turned off.


printf "\n\nTest 4: Varying the overlap \n"
for ol in $(seq $((${1:-$N}/2)))
do
   printf "Overlap: $ol out of ${1:-$N}\n"
   mpiexec -n 2 ../test1 -t1_N ${1:-$N} -snes_converged_reason -snes_type nasm -da_overlap $ol -log_view | grep 'Nonlinear solve\|Time (sec):'
done
## This is a test with 2 subdomains. So N/2 is full overlap.
## You can change the overlap with -da_overlap
## Full overlap leads to a single schwarz iteration.


printf "\n\nTest 5: Speed up\n"
M=60
for nd in 1 2 4 6 8
do
   printf "Running with n = $nd subdomains\n"
   mpiexec -n $nd ../test1 -t1_N $M -snes_converged_reason -snes_type nasm -da_overlap 2 -log_view | grep 'Nonlinear solve\|Time (sec):'
done
## The time is taken directly from the -log_view output. It doesn't seem to be wall clock time.
## We do see a pattern of decreasing time w.r.t. increasing sub domains.