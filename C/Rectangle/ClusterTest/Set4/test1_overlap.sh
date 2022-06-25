N=200

## Tw
np=4
for op in {0.0,0.05,0.10,0.15,0.20}
do
   printf "overlap = $op\n"
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex1
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex2
   mpiexec -np $np ../../maddm -N $N -op $op -problem ex3
done


# printf "\n\n - - - - - - - N=200 - - - - - - - - - \n"
# for op in {0.0,0.05,0.10,0.15,0.20}
# do
#    printf "$overlap $op\n"
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex1
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex2
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex3
# done

# ## Four Subdomains
# np=4
# printf " test with np = $np\n"
# printf " - - - - - - - N=100 - - - - - - - - - \n"
# for op in {0.0,0.05,0.10,0.15,0.20}
# do
#    printf "$overlap $op\n"
#    mpiexec -np $np ../../maddm -N $N0 -op $op -problem ex1
#    mpiexec -np $np ../../maddm -N $N0 -op $op -problem ex2
#    mpiexec -np $np ../../maddm -N $N0 -op $op -problem ex3
# done


# printf "\n\n - - - - - - - N=200 - - - - - - - - - \n"
# for op in {0.0,0.05,0.10,0.15,0.20}
# do
#    printf "$overlap $op\n"
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex1
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex2
#    mpiexec -np $np ../../maddm -N $N1 -op $op -problem ex3
# done