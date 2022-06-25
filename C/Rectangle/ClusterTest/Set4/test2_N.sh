np=4
for N in {100,200,300,400}
do
   printf "N = $N\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1
   mpiexec -np $np ../../maddm -N $N -problem ex2
   mpiexec -np $np ../../maddm -N $N -problem ex3
done