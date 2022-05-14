
snes_types="newtonls newtontr vinewtonssls fas nasm aspin"
ksp_types="cg groppcg pipecg pipecgrr pipelcg pipeprcg pipecg2 nash stcg gltr richardson chebyshev gmres tcqmr fcg pipefcg bcgs ibcgs qmrcgs fbcgs pipebcgs fbcgsr bcgsl cgs tfqmr cr pipecr lsqr preonly bicg fgmres pipefgmres minres symmlq lgmres lcd gcr pipegcr pgmres dgmres cgls"
pc_types="none jacobi pbjacobi bjacobi sor lu mg eisenstat ilu icc cholesky asm gasm ksp redundant mat cp redistribute svd gamg kaczmarz hmg lmvm"

N=20


for k in $ksp_types
do
   for p in $pc_types
   do
      printf "\nksp: $k, pc: $p\n"
      ../test1 -t1_N ${1:-$N} -snes_type newtonls -ksp_type $k -pc_type $p -snes_converged_reason -log_view | grep 'Trying\|Nonlinear\|error\|Time (sec):' | awk '{print $NF}'
   done
   printf " = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = \n"
done