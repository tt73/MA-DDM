## These tests are for the NASM with a different nonlinear solver.
N=100
np=2

# printf "Newtonls with high-tolerance settings: \n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_linesearch_type bt  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
# mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_linesearch_type bt  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
# mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_linesearch_type bt  -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Single Iteration Newton w/ ksp at 0.001 tol.

# printf "\nRichardson:\n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type nrichardson -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
## Richardson is the most basic. It just takes a step in the direction -F(u^k). It doesn't converge.

printf "\nNGMRES with SIN NPC:\n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_rtol 1e-1 -sub_npc_snes_type newtonls -sub_npc_snes_max_it 1 -sub_npc_ksp_rtol 1e-2 -snes_monitor
# mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type ngmres -sub_snes_rtol 1e-1 -sub_npc_snes_type newtonls -sub_npc_snes_max_it 1 -sub_npc_ksp_rtol 1e-2 -snes_monitor
# mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type ngmres -sub_snes_rtol 1e-1 -sub_npc_snes_type newtonls -sub_npc_snes_max_it 1 -sub_npc_ksp_rtol 1e-2 -snes_monitor
## NGMRES tries to miminize the next step using a linear combination of the previous steps.
## It chooses the directions based on -F(u^k) which is awful.
## It needs a right preconditioner to guide the direction selection. Trust region works as a NPC.
## It's not good as an inner solver.

printf "\nTR:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type newtontr -sub_snes_trtol 1e-12 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_type newtontr -sub_snes_trtol 1e-12 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 -snes_monitor
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_type newtontr -sub_snes_trtol 1e-12 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1 -snes_monitor
## Its slow but it converges without NPC
## This has 8 parameters so its annoying to tune.

exit 0

printf "\nFAS:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type fas -sub_fas_coarse_snes_linesearch_type l2 -sub_snes_max_it 1 -sub_fas_coarse_pc_type eisenstat format
## Converges without NPC

printf "\nAnderson with LS NPC:\n"
lsnpc=' -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1'
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type anderson -sub_snes_max_it 1 $lsnpc format

printf "\nAnderson with FAS NPC:\n"
fasnpc=' -sub_npc_snes_type fas -sub_npc_fas_coarse_snes_linesearch_type l2 -sub_npc_fas_coarse_pc_type eisenstat -sub_npc_fas_coarse_snes_max_it 1'
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type anderson -sub_snes_max_it 1 $fasnpc format

printf "\nNGMRES with LS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_max_it 1 $lsnpc format
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_npc_snes_type newtonls -sub_npc_pc_type eisenstat -sub_npc_snes_linesearch_type l2 -sub_npc_snes_max_it 1 -sub_snes_max_it 1 -snes_converged_reason -snes_monitor -snes_view
## NGMRES doensn't converge without a preconditioner
## This is good

printf "\nNGMRES with FAS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_max_it 1 $fasnpc format
## This is good

# printf "\nNGMRES with TR NPC:\n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type ngmres -sub_snes_max_it 1 -sub_npc_snes_type newtontr -sub_npc_pc_type eisenstat -sub_npc_snes_max_it 1 format
## This is bad. Trustregion is not a good preconditioner.


printf "\nQN with LS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type qn -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 1 $lsnpc format
## It's not fast, but it has fewer DDM iterations.

printf "\nQN with FAS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_type qn -sub_snes_linesearch_type l2 -sub_snes_linesearch_max_it 1 $fasnpc format
## It's not fast, but it has fewer DDM iterations.

