for N in 100 150 200 250 300 350 400 450 500
do
   printf "N = $N:\n"
   ../../maddm -N $N -problem ex1
   ../../maddm -N $N -problem ex2
   ../../maddm -N $N -problem ex3
done
