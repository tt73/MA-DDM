#!/bin/bash

# N = Number of interior points  
N=20
# Print the Newton iterations with -snes_monitor
# Print the SNES settings - tolerances, preconditioner, newton type, etc. with -snes_view
echo "Running code serially with $N^2 points on a square:"
../test1 -t1_N ${1:-$N} -snes_monitor -snes_view



# Set specific x & y discrete points with -da_grid_x or -da_grid_y
echo "\n\nRunning code on rectangular domain [-.5, 5] x [-1, 1]: "
../test1 -da_grid_x 15 -da_grid_y 60 -t1_Lx 0.5


echo "\n\nRunning code with ASPIN: "
../test1 -t1_N ${1:-$N} -snes_monitor -snes_view -snes_type aspin
