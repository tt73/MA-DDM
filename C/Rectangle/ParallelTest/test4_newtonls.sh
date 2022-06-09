N=100
np=4
ol=10

printf "base of comparison - GMRES + ILU:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## Tests 1, 2, and 3 were run with this type of subsolve.
## It's performs well. Almost works for all N and Nd.

## PC TYPES ============================
## none jacobi bjacobi sor lu mg eisenstat ilu asm gasm gamg
## I removed ones that don't work with GMRES from the list, like icc.


# printf "\nGMRES + multigrid:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type mg -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## Not that good.

# printf "\nGMRES + jacobi:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type jacobi -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## This one works buts its a bit slower.

# printf "\nGMRES + SOR:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type sor -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## This one is decent.

# printf "\nGMRES + eisenstat:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## This works but its actually just the same as ILU.


## KSP TYPES ============================
## richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly qcg bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls
##

# printf "\nPFGMRES + ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type pipefgmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## pipelined flexible GMRES
# ## its close to GMRES speed

# printf "\nRichardson + ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type richardson -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ## It works but it's actually quite a bit slower.

# printf "\nChebyshev + ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type chebyshev -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ##

# printf "\nGCR+ ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gcr -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# ##

# printf "\nDGMRES + ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type dgmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

# printf "\nLGMRES + ILU:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type lgmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## Linesearch types:
## shell basic l2 bt nleqerr cp ncglinear

printf "\nlinesearch bt with FAS NPC:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type bt -sub_npc_snes_type fas  -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
#

printf "\nlinesearch l2:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_npc_snes_type fas  -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
# I tried changing the linesearch method to a lower order method.
# Number of iterations doesn not change.
# It's actually a bit faster than bt used in the base settings in terms of wtime.

printf "\nlinesearch basic:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type basic -sub_npc_snes_type newtontr -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "\nlinesearch nleqerr:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type nleqerr -sub_npc_snes_type fas -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "\nlinesearch ncglinear:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type ncglinear -sub_npc_snes_type fas -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'


# printf "\ntolerance 4e-4 (h^2):\n"
# tol=4e-4
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## By default, the SNES stops when the relative residual drops below 1e-8.
## I can up the tolerance to h^2 and it will converge faster (iterations fell from 30 to 9).
## However, the final error also larger (inf error increases from 2.5e-3 to 1.3e-2).

# printf "\nAnother tolerance test 1e-6:\n"
# tol=1e-6
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## You can trade some iterations for accuracy.
## The number of iterations decreases from 30 to 21
## The error increases from 2.48e-3 to 2.50e-3.

# printf "\nStricter tolerance 1e-10:\n"
# tol=1e-10
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## I can change the tolerance to 1e-10 which is even more strict than the default.
## The number of iterations increases and the accuracy slightly improves.


# printf "\nCapped linesearch iterations: \n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 5 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

# printf "\nCapped sub_snes iterations: \n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 5 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## We keep GMRES with ILU.
## Switch BT linesearch to L2.
## Cap the sub_snes iteration with a low number for the early iterations.

# printf "\nCapped sub_snes to only one iteration: \n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This looks like the fastest way to do newton linesearch on the subdomain.

