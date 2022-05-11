#!/bin/bash
#
# Run this code with `bash basic_tests.sh`
# You can pass an integer as an arg to change N. 

# N = Number of interior points
N=20
printf "Running code serially with $N^2 points on a square:\n"
../test1 -t1_N ${1:-$N} -snes_monitor -snes_type newtonls
# Print the Newton iterations with -snes_monitor


# Set specific x & y discrete points with -da_grid_x or -da_grid_y
printf "\n\nRunning code on rectangular domain [-.5, 5] x [-1, 1]:\n"
../test1 -t1_Lx 0.5
# This works just fine


# Convergence test
# Run the program for N = 10, 20, ...
printf "\n\nConvergence test:\n"
for i in {10..50..10}
do
   # ../test1 -t1_N $i -snes_type aspin
done
# the error descreases as expected
# the epsilon and width dynamically with N


#
printf "\n\nShow list of availables SNES types:\n"
../test1 -help | grep snes_type

snes_types="newtonls newtontr nrichardson ksponly ksptransposeonly vinewtonssls ngmres qn ncg fas nasm anderson aspin"
for s in $snes_types
do
   printf "\nTrying out $s:\n"
   ../test1 -t1_N ${1:-$N} -snes_type $s
done
# Converge: newtonls, newtontr, vinewtonssls, fas, nasm, aspin
# Diverge: nrichardson, ksponly, ksptransposeonly, ngmres, qn, ncg, anderson


printf "\n\nTiming the code:\n"
../test1 -t1_N 50  -log_view | grep 'error\|problem\|Time (sec):'