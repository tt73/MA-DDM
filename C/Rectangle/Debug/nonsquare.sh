## This is a unit test for when Nx != Ny.

## We know square works. This is h = 0.5.
printf "This is a working example with Nx = Ny = 3. Notice the error: \n"
../maddm -ymin -2 -ymax 2 -xmin -1 -xmax 1 -Nx 3 -Ny 3

## Try running code on [-1,1] x [-2,2] with hx = hy, ie. Ny = 2*Nx+1.
## Small example with hx = hy = 0.5.
## Stencil width = 1 actually works
Nx=3
Ny=7
printf "\nRunning on small problem Nx = $Nx, Ny = $Ny with width=1... it should work\n"
../maddm -ymin -2 -ymax 2 -xmin -1 -xmax 1 -Nx $Nx -Ny $Ny -width 1

printf "\nRunning on Nx = $Nx, Ny = $Ny now with width=2... \n"
../maddm -ymin -2 -ymax 2 -xmin -1 -xmax 1 -Nx $Nx -Ny $Ny -width 2
printf "If you see *Iters: 9999 then it failed.\n"

## Stencil width > 1 seems to fail.
## ComputeWeights might possibly be the issue. I tried 2nd order but it still doesn't work.
## The most likely culprits are ComputeSDD or ComputeProjectionIndeces.
## Can I isolate the issue somehow?

printf "\n Debugging \n\n "
../maddm -ymin -2 -ymax 2 -xmin -1 -xmax 1 -Nx $Nx -Ny $Ny -width 2 -debug -snes_max_it 1
