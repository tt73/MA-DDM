N=100
np=4
ol=10

## base
printf "base of comparison:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

# printf "multigrid pc:\n"
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type mg -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "linesearch l2:\n"
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "tolerance 4e-4 (h^2):\n"
tol=4e-4
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "tolerance 1e-6:\n"
tol=1e-6
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "tolerance 1e-7:\n"
tol=1e-7
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "tolerance 1e-10:\n"
tol=1e-10
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

printf "tolerance npc newtonls:\n"
tol=1e-10
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -npc_snes_type newtonls -sub_ksp_type gmres -sub_pc_type ilu -snes_rtol $tol -snes_atol $tol -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

