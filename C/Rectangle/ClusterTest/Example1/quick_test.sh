#!/bin/bash -l


# mpiexec -np 4 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-4 -sub_ksp_rtol 1e-3 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-8

# mpiexec -np 4 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-3 -sub_ksp_rtol 1e-3 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-8

# mpiexec -np 4 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-4 -sub_ksp_rtol 1e-4 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-7

# mpiexec -np 4 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-6 -sub_ksp_rtol 1e-6 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-6 -sol



mpiexec -np 6 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-2 -sub_ksp_rtol 1e-2 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-8

mpiexec -np 6 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-3 -sub_ksp_rtol 1e-3 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-8

mpiexec -np 6 ../../maddm -N 200 -problem ex1 -sub_snes_rtol 1e-4 -sub_ksp_rtol 1e-4 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-8