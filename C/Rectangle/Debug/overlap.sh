
## In theory, the DDM convergence rate should increase with respect to domain overlap.
## The iteration number should decrease as `da_overlap` (an integer amount of overlap in all directions) increases.

# load petsc
module load gnu8 mpich petsc/3.12.0

clear

## Settings
N=32     ## N=2^5
prob=ex1 ## gaussian [-1,1]^2
h=$(echo "scale = 8; 2/($N+1)" | bc)
h2=$(echo "scale = 8; $h*$h" | bc)
np=4     ## four subdomains


file=op.out
rm -f $file

for op in {0..32}
do
   printf "op = $op\n" >> $file
   mpiexec -np $np ../maddm -N $N -problem ex1 -sin  -snes_converged_reason -snes_atol $h -snes_rtol 1e-99 -sub_ksp_type preonly -sub_pc_type lu -sub_snes_linesearch_type basic -da_overlap $op >> $file
done

