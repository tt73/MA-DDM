N=300

for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
   mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
   mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
done


for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
   mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
   mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
done
