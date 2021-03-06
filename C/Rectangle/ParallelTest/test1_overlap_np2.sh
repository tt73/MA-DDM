N0=100
N1=200

## Two subdomains
np=2

printf " - - - - - - - N=100 - - - - - - - - - \n"
## 0% overlap
printf "0%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 0 -t1_problem ex12

## 5% overlap
printf "5%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 2 -t1_problem ex12

## 10% ovrelap
printf "10%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 5 -t1_problem ex12

## 15% ovrelap
printf "15%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 7 -t1_problem ex12

## 20% ovrelap
printf "20%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N0 -da_overlap 10 -t1_problem ex12



printf "\n\n - - - - - - - N=200 - - - - - - - - - \n"
printf "0%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 0 -t1_problem ex12

## 5% overlap
printf "5%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 5 -t1_problem ex12

## 10% ovrelap
printf "10%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 10 -t1_problem ex12

## 15% ovrelap
printf "15%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 15 -t1_problem ex12

## 20% ovrelap
printf "20%% overlap\n"
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex10
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex11
mpiexec -np $np ../test1 -t1_N $N1 -da_overlap 20 -t1_problem ex12
