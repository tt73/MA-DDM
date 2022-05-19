#!/bin/bash

# We let dmda manage the domain decomposition. Petsc will automatically handle the decomposition.
# Sometimes, the petsc code returns an error if the MPI size is not compatible with grid size.
# We need to figure what's acutally happening with the domain decomposition.
printf "Test 1: DM View\n"
printf " --------------------------------------------\n"
mpiexec -n 2 ../test1 -t1_N 10 -dm_view -snes_view
printf "   This is output for NASM with N = 10, and np = 2\n"
printf "   It appears that the default overlap is 1.\n"
printf "   The domain is split in top and bottom halves.\n"
printf " --------------------------------------------\n"
mpiexec -n 3 ../test1 -t1_N 10 -dm_view -da_overlap 2
printf "   This is output for NASM with N = 10, and np = 3\n"
printf "   There is no overlap. The domain is split in top and bottom halves.\n"
printf " --------------------------------------------\n"

