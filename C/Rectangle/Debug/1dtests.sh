## This set of tests is to make sure the 1D code is working


## Just to review the options for maddm
printf "Printing options for maddm:\n"
../maddm -help | grep maddm
printf "\n\n"

## To use the 1D version, you add the option -dim 1.
## Check the plot to make sure its correct.
printf "Running 1D code with N = 100\n"
../maddm -dim 1 -N 100 -sol
python3 plot1d.py -n 100 # save solution plot
## The error should be around 1e-5

## Problem 2 should also have low error
printf "\n\nNow doing problem 2\n"
../maddm -dim 1 -N 100 -problem ex2
## 1e-4


## Problem 3
printf "\n\nNow doing problem 3\n"
../maddm -dim 1 -N 100 -problem ex3
## 1e-6


## Problem 4
printf "\n\nNow doing problem 4\n"
../maddm -dim 1 -N 100 -problem ex4