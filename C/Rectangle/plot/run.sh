module load python
module load gnu8 mpich petsc/3.12.0

N=300

../maddm -N $N -problem ex1 -sol
mv load_u.m load_u1.m
mv load_exact.m load_exact1.m

../maddm -N $N -problem ex2 -sol
mv load_u.m load_u2.m
mv load_exact.m load_exact2.m

../maddm -N $N -problem ex3 -sol
mv load_u.m load_u3.m
mv load_exact.m load_exact3.m

../maddm -N $N -problem ex4 -sol
mv load_u.m load_u4.m
mv load_exact.m load_exact4.m

