
N=301

d01='-t1_xmin 0 -t1_ymin 0'

printf "Domain test with N = $N\n"

np=1
printf "\nOne domain\n"
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

np=2
ox=0
oy=15
printf "\nTwo-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

np=4
ox=15
oy=15
printf "\nFour-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

np=6
ox=15
oy=10
printf "\nSix-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1

np=9
ox=10
oy=10
printf "\nNine-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 $d01 -sub_ksp_type gmres -sub_pc_type eisenstat -sub_snes_linesearch_type l2 -sub_snes_max_it 1