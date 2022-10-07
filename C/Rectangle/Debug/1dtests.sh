## This set of tests is to make sure the 1D code is working


## Just to review the options for maddm
printf "Printing options for maddm:\n"
../maddm -help | grep maddm
printf "\n\n"

## To use the 1D version, you add the option -dim 1.
## Check the plot to make sure its correct.
printf "Running 1D code with N = 100\n"
../maddm -dim 1 -N 100 -sol -width 1
python3 plot1d.py -n 100 # save solution plot
## The error should be around 1e-5

## Problem 2 should also have low error
printf "\n\nNow doing problem 2\n"
../maddm -dim 1 -N 100 -problem ex2 -width 1
## 1e-4


## Problem 3
printf "\n\nNow doing problem 3\n"
../maddm -dim 1 -N 100 -problem ex3 -width 1
## 1e-6


## Problem 4
printf "\n\nNow doing problem 4\n"
../maddm -dim 1 -N 100 -problem ex4 -width 1


printf "\n\nRepeat problem 1 with d = 2 with FD Jacobian\n"
../maddm -dim 1 -N 100 -sol -width 2 -sub_snes_fd
## this should work
## if it doesn't, there's something wrong with F

printf "\n\nRepeat problem 1 with d = 2 with handcoded Jacobian\n"
../maddm -dim 1 -N 100 -sol -width 2
## its way worse than FD for some odd reason
## but the error looks fine

printf "\n\nRepeat problem 1 with d = 3 with handcoded Jacobian\n"
../maddm -dim 1 -N 100 -sol -width 3

printf "\n\nRepeat problem 1 with d = 3 with FD Jacobian\n"
../maddm -dim 1 -N 100 -sol -width 3 -sub_snes_fd


