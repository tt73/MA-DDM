## Problem 3 aka Example 12



N=300
ol=15
np=4
printf "\nN = $N with ol = $ol\n"

## One-step Newton fails to converge
# mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex12

## Mixed NASM will converge. Takes around 4 min, 59 iterations.
## Its doing one-step Newton on the 3 regular domains and doing 3-step cubic bt on the last one.
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 -t1_mixed -snes_monitor

## NKS takes around 2 min, 58 iterations.
mpiexec -np $np ../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 -snes_type newtonls -snes_monitor