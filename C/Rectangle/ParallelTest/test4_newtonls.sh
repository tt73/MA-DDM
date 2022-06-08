N=100
np=4
ol=10

printf "base of comparison:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## Tests 1, 2, and 3 were run with this type of subsolve.
## It's performs well. Almost works for all N and Nd.

# printf "PC - multigrid:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type mg -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This one stalls and doesn't converege.

printf "\nPC - jacobi:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type jacobi -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This one works buts its a bit slower.

printf "\nPC - block jacobi:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type bjacobi -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This one works and its a bit faster than ilu sometimes.

printf "\nlinesearch l2:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## I tried changing the linesearch method to a lower order method.
## Number of iterations doesn not change.
## It's actually a bit faster than bt used in the base settings in terms of wtime.

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

# printf "\nKSP - CG:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type cg -sub_pc_type sor -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## Conjugate gradient is not working.

# printf "\nKSP - Richardson:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type richardson -sub_pc_type gasm -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## It works but it's actually quite a bit slower.

# printf "\nKSP - genralized conjugate residual:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type richardson -sub_pc_type gasm -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This method is for non-symmetric matrices.
## Too slow.

# printf "\nKSP - deflated gmres:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type dgmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This method is a variation of gmres.
## Its


printf "\nCapped linesearch iterations: \n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 5 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "\nCapped sub_snes iterations: \n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 5 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## We keep GMRES with ILU.
## Switch BT linesearch to L2.
## Cap the sub_snes iteration with a low number for the early iterations.

printf "\nCapped sub_snes to only one iteration: \n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
## This looks like the fastest way to do newton linesearch on the subdomain.

