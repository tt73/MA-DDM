

N0=50
N1=200
# N2=200

## Two subdomains
np=2

## 0% overlap
# mpiexec -np $np ../test1 -t1_N $N1 -snes_converged_reason -da_overlap 0

## 5% overlap
# mpiexec -np $np ../test1 -t1_N $N1 -snes_converged_reason  -da_overlap 5

## 10% ovrelap
# mpiexec -np $np ../test1 -t1_N $N1 -snes_converged_reason -da_overlap 10

## 15% ovrelap
# mpiexec -np $np ../test1 -t1_N $N1 -snes_converged_reason -da_overlap 15

## 20% ovrelap
# mpiexec -np $np ../test1 -t1_N $N1 -snes_monitor -snes_converged_reason -snes_view -log_view -da_overlap 20


# mpiexec -np $np ../test1 -t1_N $N0 -snes_monitor -snes_converged_reason -snes_view -log_view -da_overlap 0

# mpiexec -np $np ../test1 -t1_N $N0 -snes_monitor -snes_converged_reason -snes_view -log_view -da_overlap 10  -npc_snes_type newtonls

mpiexec -np $np ../test1 -t1_N $N1 -snes_monitor -snes_converged_reason -log_view -snes_view -da_overlap 5 -sub_ksp_type gmres -sub_pc_type ilu