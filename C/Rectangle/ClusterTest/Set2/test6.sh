N=200
np=4

printf "\nNASM:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex2 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1
mpiexec -np $np ../../maddm -N $N -problem ex3 -sub_snes_rtol 1e-1 -sub_ksp_rtol 1e-1

printf "\nLS:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type newtonls -snes_linesearch_type l2 -ksp_type pipefgmres -ksp_rtol 1e-2

printf "\nASPIN:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type aspin -npc_sub_pc_type ilu
mpiexec -np $np ../../maddm -N $N -problem ex2 -snes_type aspin -npc_sub_pc_type ilu
mpiexec -np $np ../../maddm -N $N -problem ex3 -snes_type aspin -npc_sub_pc_type ilu


printf "\nFAS:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type fas -fas_coarse_snes_linesearch_type l2

# printf "\nNGMRES:\n"
# mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres -snes_ngmres_select_type difference




printf "\nNGMRES with Newton LS NPC:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres -snes_npc_side right -npc_snes_type newtonls -npc_snes_linesearch_type basic -npc_snes_max_it 1


printf "\nNGMRES with FAS NPC:\n"
fasnpc='-npc_snes_type fas -npc_fas_coarse_snes_linesearch_type l2 -npc_snes_max_it 1 -npc_fas_levels_snes_type newtonls -npc_fas_levels_snes_max_it 6'
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type ngmres $fasnpc


printf "\nFAS * NGMRES composite:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type composite -snes_composite_type multiplicative -snes_composite_sneses fas,ngmres

printf "\nFAS + Newton composite:\n"
mpiexec -np $np ../../maddm -N $N -problem ex1 -snes_type composite -snes_composite_type additiveoptimal -snes_composite_sneses newtontr,newtonls