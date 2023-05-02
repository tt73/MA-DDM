## I want to see what the initial guesses look like.
## I also want to see what they look like after one iteration of Newton.


N = 32

../../maddm -N $N -snes_max_it 1 -initial zeros -debug -sol

# ../maddm -N $N -snes_max_it 1 -initial random -debug -sol