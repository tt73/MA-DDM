N=300
printf "N = $N\n"
for np in {1,9}
# for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex2 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex3
done

N=400
printf "N = $N\n"
for np in {1,9}
# for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex2 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex3
done

N=500
printf "N = $N\n"
for np in {1,9}
# for np in {1,2,4,6,8,9}
do
   printf "Nd = $np\n"
   mpiexec -np $np ../../maddm -N $N -problem ex1 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex2 -fast
   mpiexec -np $np ../../maddm -N $N -problem ex3
done