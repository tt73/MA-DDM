N=300
for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1
   mpiexec -np $np ../../maddm -N $N -problem ex2
   mpiexec -np $np ../../maddm -N $N -problem ex3
done