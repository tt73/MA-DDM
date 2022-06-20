N0=100
N1=200

## Four subdomains
np=4

printf " test with np = $np\n"
printf " - - - - - - - N=100 - - - - - - - - - \n"
## 0% overlap
printf "0%% overlap\n"
mpiexec -np $np ../../maddm -N $N0 -op 0.0 -problem ex1
mpiexec -np $np ../../maddm -N $N0 -op 0.0 -problem ex2
mpiexec -np $np ../../maddm -N $N0 -op 0.0 -problem ex3

## 5% overlap
printf "5%% overlap\n"
mpiexec -np $np ../../maddm -N $N0 -op 0.05 -problem ex1
mpiexec -np $np ../../maddm -N $N0 -op 0.05 -problem ex2
mpiexec -np $np ../../maddm -N $N0 -op 0.05 -problem ex3

## 10% ovrelap
printf "10%% overlap\n"
mpiexec -np $np ../../maddm -N $N0 -op 0.1 -problem ex1
mpiexec -np $np ../../maddm -N $N0 -op 0.1 -problem ex2
mpiexec -np $np ../../maddm -N $N0 -op 0.1 -problem ex3

## 15% ovrelap
printf "15%% overlap\n"
mpiexec -np $np ../../maddm -N $N0 -op 0.15 -problem ex1
mpiexec -np $np ../../maddm -N $N0 -op 0.15 -problem ex2
mpiexec -np $np ../../maddm -N $N0 -op 0.15 -problem ex3

## 20% ovrelap
printf "20%% overlap\n"
mpiexec -np $np ../../maddm -N $N0 -op 0.2 -problem ex1
mpiexec -np $np ../../maddm -N $N0 -op 0.2 -problem ex2
mpiexec -np $np ../../maddm -N $N0 -op 0.2 -problem ex3



printf "\n\n - - - - - - - N=200 - - - - - - - - - \n"
printf "0%% overlap\n"
mpiexec -np $np ../../maddm -N $N1 -op 0.0 -problem ex1
mpiexec -np $np ../../maddm -N $N1 -op 0.0 -problem ex2
mpiexec -np $np ../../maddm -N $N1 -op 0.0 -problem ex3

## 5% overlap
printf "5%% overlap\n"
mpiexec -np $np ../../maddm -N $N1 -op 0.05 -problem ex1
mpiexec -np $np ../../maddm -N $N1 -op 0.05 -problem ex2
mpiexec -np $np ../../maddm -N $N1 -op 0.05 -problem ex3

## 10% ovrelap
printf "10%% overlap\n"
mpiexec -np $np ../../maddm -N $N1 -op 0.10 -problem ex1
mpiexec -np $np ../../maddm -N $N1 -op 0.10 -problem ex2
mpiexec -np $np ../../maddm -N $N1 -op 0.10 -problem ex3

## 15% ovrelap
printf "15%% overlap\n"
mpiexec -np $np ../../maddm -N $N1 -op 0.15 -problem ex1
mpiexec -np $np ../../maddm -N $N1 -op 0.15 -problem ex2
mpiexec -np $np ../../maddm -N $N1 -op 0.15 -problem ex3

## 20% ovrelap
printf "20%% overlap\n"
mpiexec -np $np ../../maddm -N $N1 -op 0.20 -problem ex1
mpiexec -np $np ../../maddm -N $N1 -op 0.20 -problem ex2
mpiexec -np $np ../../maddm -N $N1 -op 0.20 -problem ex3
