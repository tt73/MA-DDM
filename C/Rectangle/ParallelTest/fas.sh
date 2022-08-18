## FAS is short for Full approximation multigrid scheme and it has many customization options.
## Run `./maddm -snes_type fas -help | grep fas` to see all of the runtime options.
## This is the manual page: https://petsc.org/main/docs/manualpages/SNESFAS/SNESFAS/
N=100 # grid is N by N
np=4 # number of MPI processes



printf "Test 1: Using FAS directly\n"
printf "default serial: "
../maddm -N $N -snes_type fas -snes_converged_reason
## You can set the main snes_type to FAS and run it.
printf "additive serial: "
../maddm -N $N -snes_type fas -snes_converged_reason -snes_fas_type additive
## You can change multiplicative to additive for parallelization
printf "additive more cycles serial: "
../maddm -N $N -snes_type fas -snes_converged_reason -snes_fas_type additive -snes_fas_cycles 3
## You can add more cycles. Not sure what this does.
printf "additive parallel: "
mpiexec -np $np ../maddm -N $N -snes_type fas -snes_converged_reason -snes_fas_type additive
## You can run fas in parallel. Even though its additive, the parallelization comes from using a block-Jacobi precenditioner on the coarse level.
printf "limit coarse solve iterations: "
mpiexec -np $np ../maddm -N $N -snes_type fas -snes_converged_reason -fas_coarse_snes_max_it 5
## You can run fas in parallel. Even though its additive, the parallelization comes from using a block-Jacobi precenditioner on the coarse level.
printf "two-level FAS in parallel: "
mpiexec -np $np ../maddm -N $N -snes_type fas -snes_fas_levels 2 -fas_levels_snes_type newtonls -fas_levels_snes_max_it 5 -snes_converged_reason -fas_coarse_snes_max_it 5
## You can run fas in parallel. Even though its additive, the parallelization comes from using a block-Jacobi precenditioner on the coarse level.



## FAS can be used as a nonlinear precenditioner.
printf "\n\nTest 2: FAS as an NPC\n"
printf "serial nrichardson with fas npc: "
../maddm -N $N -snes_type nrichardson -npc_snes_type fas -snes_converged_reason
## You can let nonlinear richardson be the main method and let FAS be the preconditioner.
printf "serial ngmres with fas npc: "
../maddm -N $N -snes_type ngmres -npc_snes_type fas -snes_converged_reason
##
printf "parallel ngmres with fas npc: "
mpiexec -np $np ../maddm -N $N -snes_type ngmres -npc_snes_type fas -snes_converged_reason
printf "parallel ngmres with fas npc + more cycles: "
mpiexec -np $np ../maddm -N $N -snes_type ngmres -npc_snes_type fas -npc_snes_fas_cycles 3 -npc_snes_fas_smoothup 2 -snes_converged_reason


printf "\n\nTest 3: NASM with FAS on the subdomains\n"
printf "NASM FAS with defualt settings: "
mpiexec -np $np ../maddm -N $N -snes_type nasm -da_overlap 5 -sub_snes_type fas -snes_converged_reason -snes_rtol 1e-4
## NASM FAS with default settings is really slow. The last iterations take a long time to run.
## We can cap the residue tolerance to 1e-4 which gives an accurate answer.
printf "Limited inner iterations: "
mpiexec -np $np ../maddm -N $N -snes_type nasm -da_overlap 5 -sub_snes_type fas -sub_snes_max_it 1 -sub_fas_coarse_snes_max_it 5 -snes_converged_reason -snes_rtol 1e-4

printf "NASM with NGMRES-FAS on the subdomains: "
mpiexec -np $np ../maddm -N $N -snes_type nasm -da_overlap 5 -sub_snes_type ngmres -sub_npc_snes_type fas -snes_converged_reason -snes_rtol 1e-4


printf "NASM with NGMRES-FAS with limited inner iterations: "
mpiexec -np $np ../maddm -N $N -snes_nasm_damping 0.9 -da_overlap 5 -sub_snes_type ngmres -sub_npc_snes_type fas -sub_snes_max_it 1 -sub_npc_fas_coarse_snes_max_it 5 -snes_converged_reason