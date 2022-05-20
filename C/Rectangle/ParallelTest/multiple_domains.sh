#!/bin/bash

# We let dmda manage the domain decomposition. Petsc will automatically handle the decomposition.
# Sometimes, the petsc code returns an error if the MPI size is not compatible with grid size.
# We need to figure what's acutally happening with the domain decomposition.

printf "Test 1: Ghost points\n"
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
##
## The local indeces that we get from DMLocalInfo don't show the effect of the overlap.
## The ghost points are fully determined by the stencil width.

printf "Test 2: Odd splitting\n"
printf " --------------------------------------------\n"
mpiexec -n 3 ../test1 -t1_N 10 -t1_debug -da_overlap 1 -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 3\n"
printf "* The domain is split into a vertical stack of 3 subdomains\n"
printf " --------------------------------------------\n"
mpiexec -n 5 ../test1 -t1_N 10 -da_overlap 1 -t1_debug -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 5\n"
printf "* Similar to before, the subdomains are vertical stacks.\n"
printf " --------------------------------------------\n"
mpiexec -n 7 ../test1 -t1_N 10 -da_overlap 1 -t1_debug -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 7\n"
printf "* Similar to before, the subdomains are vertical staks.\n"
printf " --------------------------------------------\n"
mpiexec -n 9 ../test1 -t1_N 10 -da_overlap 1 -t1_debug -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 9\n"
printf "* Finally, we get a tiled 3 by 3 splitting.\n"
printf "* This takes the longest to compute in wallclock time but it takes fewer iterations to converge.\n"
printf " --------------------------------------------\n\n"
##
## Odd number of subdomains usually results in vertical stacks.
## np = 9 is an exception. This results in a 3 by 3 tiled configuration.


printf "Test 3: Even splitting\n"
printf " --------------------------------------------\n"
mpiexec -n 4 ../test1 -t1_N 10 -t1_debug -da_overlap 1 -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 4\n"
printf "* The domain is split into 2 by 2.\n"
printf " --------------------------------------------\n"
mpiexec -n 6 ../test1 -t1_N 10 -da_overlap 1 -t1_debug -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 6\n"
printf "* Domain is split into 3 by 2.\n"
printf " --------------------------------------------\n"
mpiexec -n 12 ../test1 -t1_N 10 -da_overlap 1 -t1_debug -t1_width 1 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 12\n"
printf "* Domain is split into 4 by 3.\n"
printf " --------------------------------------------\n\n"
##
## For even number of subdomains, PETSc tries to make a tiling configuration as close to a square as possible.
## This is apparent from the example with 12 subdomains which got split into 3 stacks of 4.
## PETSc opts to make number of rows larger than number of columns when it comes to even subdomains.


printf "Test 4: Errors\n"
printf " --------------------------------------------\n"
mpiexec -n 5 ../test1 -t1_N 10 -t1_debug -da_overlap 1 -t1_width 3 -snes_converged_reason
printf "* This is output for NASM with N = 10, and np = 5\n"
printf "* Each subdomain is 2 nodes high so a width of 3 gives an error.\n"
printf " --------------------------------------------\n\n"
##
## The local indeces that we get from DMLocalInfo don't show the effect of the overlap.
## The ghost points are fully determined by the stencil width.