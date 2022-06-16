
N=512

d01='-t1_xmin 0 -t1_ymin 0'

printf "Domain test with N = $N\n"

np=1
printf "\nOne domain\n"
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

np=2
ox=0
oy=25
printf "\nTwo-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

np=4
ox=25
oy=25
printf "\nFour-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

np=6
ox=25
oy=17
printf "\nSix-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

np=9
ox=17
oy=17
printf "\nNine-domain, ox = $ox, oy = $oy\n"
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

# np=12
# ox=17
# oy=12
# printf "\nTwelve-domain, ox = $ox, oy = $oy\n"
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

# np=15
# ox=17
# oy=10
# printf "\nFifteen-domain, ox = $ox, oy = $oy\n"
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2

# np=16
# ox=12
# oy=12
# printf "\nFifteen-domain, ox = $ox, oy = $oy\n"
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex10 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex11 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2
# mpiexec -np $np ../../test1 -t1_N $N -da_overlap_x $ox -da_overlap_y $oy -t1_problem ex12 -t1_xmin 0 -t1_ymin 0 -snes_type newtonls -snes_linesearch_type l2