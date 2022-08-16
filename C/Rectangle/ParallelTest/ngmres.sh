N=100
Nd=2

## Regular SIN
mpiexec -np $Nd ../maddm -N $N -problem ex1 -sin


## I don't think its possible to have an NPC for NASM, so NASM-NGMRES is not doable

## NGMRES with NASM NPC is possible
mpiexec -np $Nd ../maddm -N $N -problem ex1 -snes_type ngmres -npc_snes_type nasm -npc_snes_nasm_type restrict

## NGMRES with SIN NPC
mpiexec -np $Nd ../maddm -N $N -problem ex1 -snes_type ngmres -npc_snes_type nasm -npc_snes_nasm_type restrict -npc_sub_snes_max_it 1 -npc_sub_ksp_rtol 1e-1 -npc_sub_ksp_type dgmres -npc_sub_pc_type eisenstat
## this

