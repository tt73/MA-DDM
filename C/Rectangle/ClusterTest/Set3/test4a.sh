## This is the initial test to make the local solver faster.
## Want to find out how far we can relax the newton and kyrlov tolerances before the method suffers.

N=200
np=4

for snestol in {1.e-4,1.e-3,1.e-2,1.e-1}
do
   for ksptol in {1.e-4,1.e-3,1.e-2,1.e-1}
   do
      printf "snestol = $snestol, ksptol = $ksptol\n"
      mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol
      mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol
      mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_rtol $snestol -sub_ksp_rtol $ksptol
   done
done