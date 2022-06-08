N0=100
N1=200

## Two subdomains
np=2

printf " test with np = $np\n"
printf " - - - - - - - N=100 - - - - - - - - - \n"
## 0% overlap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason -da_overlap 0
printf "0%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 5% overlap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason  -da_overlap 5
printf "5%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1  -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1  -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1  -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 10% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason -da_overlap 10
printf "10%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 15% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason -da_overlap 15
printf "15%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 20% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_monitor -snes_converged_reason -snes_view -log_view -da_overlap 20
printf "20%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'



printf "\n\n - - - - - - - N=200 - - - - - - - - - \n"
printf "0%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 5% overlap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason  -da_overlap 5
printf "5%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 10% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason -da_overlap 10
printf "10%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 15% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_converged_reason -da_overlap 15
printf "15%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'

## 20% ovrelap
# mpiexec -np $np ../../test1 -t1_N $N1 -snes_monitor -snes_converged_reason -snes_view -log_view -da_overlap 20
printf "20%% overlap\n"
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex11 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
mpiexec -np $np ../../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex12 -sub_ksp_type gmres -sub_pc_type ilu -sub_snes_linesearch_type l2 -sub_snes_max_it 1 -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
