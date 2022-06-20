#!/bin/bash
# The purpose of this script is to show how the MA-DDM program can be run.
#
# Run this code with `bash basic_usage.sh`
# You can pass an integer as an arg to change N.


# N = Number of interior points
N=20

printf "Running code serially with $N^2 points on a square.\nNewton iterations: \n"
../maddm -N ${1:-$N} -snes_monitor -snes_type newtonls -sol
## Change size of problem with -N
## Change the outer nonlinear solver with -snes_type.
## Turn on option to print residues with -snes_monitor.
## Turn on option to generate matlab solution scripts with -sol.


printf "Running code on 40 by 30 grid.\n"
../maddm -Nx 30 -Ny 40 -width 1
## The size in the x and y directions can be individually with -Nx and -Ny.


printf "Manually set the stencil width and epsilon\n"
../maddm -N 20 -width 1 -eps 0.01
## The stencil width is set to the ceiling of d=(h)^-(1/3),
## and the epsilon is set to h^2. These values can be
## manually changed. However, there is an optimal d for any fixed h.


printf "\n\nRunning code on rectangular domain [-.5, .5] x [-2, 2]:\n"
../maddm -N $N -xmin -0.5 -xmax 0.5 -ymin -2.0 -ymax 2.0
## The limits of the rectangular domain can be changed with -xmin/-xmax adn -ymin/-ymax.

printf "\n\nChange the problem :\n"
printf "example 1: \n"; ../maddm -N 10
printf "example 2: \n"; ../maddm -N 10 -problem ex2
printf "example 3: \n"; ../maddm -N 10 -problem ex3
## You can switch between 3 different example problems.
## The performance of the solver varies for each one.

printf "Change examples to 1D "
../maddm -N 100 -dim 1
## Examples 1, 2, and 3 all have a 1D version.
## A toggle for 3D also exists but it is not supported.


## Change initial guess
printf "\n\nTry out initial different guesses for the Newton iteration :\n"
printf "zeros:   "; ../maddm -N 40 -snes_converged_reason -snes_type newtontr -init_type zeros | grep 'Nonlinear'
printf "random:  "; ../maddm -N 40 -snes_converged_reason -snes_type newtontr -init_type random | grep 'Nonlinear'
printf "corner:  "; ../maddm -N 40 -snes_converged_reason -snes_type newtontr -init_type corner | grep 'Nonlinear'
printf "pyramid: "; ../maddm -N 40 -snes_converged_reason -snes_type newtontr -init_type pyramid | grep 'Nonlinear'
## Zeros is the simplest, u0 = 0.
## Random sets the initial guess to be uniform random.
## Corner sets the u0 = M, where M is the largest value of g on the 4 corners.
## Pyramid sets the initial guess to be an inverted pyramid.