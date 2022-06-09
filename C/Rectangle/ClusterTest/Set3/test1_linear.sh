snes_types="newtonls newtontr vinewtonssls fas nasm aspin"
ksp_types="cg groppcg pipecg richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr lsqr preonly bicg fgmres lgmres lcd gcr pipegcr pgmres dgmres cgls"
pc_types="none jacobi sor lu mg eisenstat ilu icc cholesky ksp redundant mat cp redistribute svd gamg kaczmarz hmg lmvm"


N=100
np=4
ol=10


for k in $ksp_types
do
   for p in $pc_types
   do
      printf "\nksp: $k, pc: $p\n"
      mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex10 -sub_ksp_type $k -sub_pc_type $p -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
      mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex11 -sub_ksp_type $k -sub_pc_type $p -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
      mpiexec -np $np ../../test1 -t1_N $N -da_overlap $ol -t1_problem ex12 -sub_ksp_type $k -sub_pc_type $p -snes_converged_reason -log_view | grep '*Problem\|*Error\|WTime\|Nonlinear solve'
   done
   printf " = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n"
done