
np=2

d01='-t1_xmin 0 -t1_ymin 0'

printf "Increasing N with fixed overlap percentage. np = $np subdomains.\n"

N=100
ol=5
printf "\nN = $N with ol = $ol\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

N=200
ol=10
printf "\nN = $N with ol = $ol\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

N=300
ol=15
printf "\nN = $N with ol = $ol\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
## for reasons unknown, problem 12 gives NaN when ran with N=300, N=301 however works
mpiexec -np $np ../../test1 -t1_N 301 -da_overlap $ol -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

N=400
ol=20
printf "\nN = $N with ol = $ol\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1