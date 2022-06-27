## Single iteration newton

N=200
np=4


for ksptol in {1.e-4,1.e-3,1.e-2,1.e-1}
do
   printf "ksptol = $ksptol\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol -sub_snes_linesearch_type l2
   mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol -sub_snes_linesearch_type l2
   mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_max_it 1 -sub_ksp_rtol $ksptol -sub_snes_linesearch_type l2
done
