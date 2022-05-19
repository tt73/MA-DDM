#!/bin/bash

# We let dmda manage the domain decomposition. Petsc will automatically handle the decomposition.
# Sometimes, the petsc code returns an error if the MPI size is not compatible with grid size.
# We need to figure what's acutally happening with the domain decomposition.
printf "Test 1: Overlap\n"
printf " --------------------------------------------\n"
mpiexec -n 2 ../test1 -t1_N 10 -t1_debug -da_overlap 1 -t1_width 1
printf "* This is output for NASM with N = 10, and np = 2\n"
printf "* The stencil width is chosen to be 1.\n"
printf "* The DDM overlap is chosen to be 1.\n"
printf "* The nodes are indexed from 0 to 9 in x & y.\n"
printf "* The domain is split into top and bottom halves as indicated by the local indeces of rank 0 and 1.\n"
printf "* Notice the ghost point indeces in parantheses.\n"
printf " --------------------------------------------\n"
mpiexec -n 2 ../test1 -t1_N 10 -da_overlap 2 -t1_debug -t1_width 1
printf "* This is output for NASM with N = 10, and np = 2\n"
printf "* The width is 1 again.\n"
printf "* The DDM overlap is changed to 2.\n"
printf "* Notice that ghost points are independent of the overlap.\n"
printf " --------------------------------------------\n\n"




