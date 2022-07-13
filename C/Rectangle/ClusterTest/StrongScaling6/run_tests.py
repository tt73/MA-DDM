import numpy as np
import subprocess

nps = np.array([1,2,4,6,8,9])

# loop to create job files
for i in range(len(nps)):
   jobfile = "job{:d}".format(nps[i])
   f = open(jobfile, "w")
   f.write("#!/bin/bash -l\n")
   f.write("#SBATCH -J j{:d}\n".format(nps[i]))
   f.write("#SBATCH -p public\n")
   f.write("#SBATCH -o out{:d}\n".format(nps[i]))
   f.write("#SBATCH -p dms-cpu\n")
   f.write("#SBATCH --mail-user tt73@njit.edu\n")
   f.write("#SBATCH -A tt73\n")
   f.write("#SBATCH -t 8:0:0\n")
   f.write("#SBATCH --mem=16G\n")
   f.write("#SBATCH --nodes {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks {:d}\n".format(nps[i]))
   # f.write("#SBATCH --exclusive\n")
   f.write("module load gnu8 mpich petsc\n")
   f.write("N=300\n")
   f.write("mpirun ../../maddm -N $N -problem ex1 -snes_type newtonls -ksp_rtol 1e-2 -snes_linesearch_order 2 -ksp_type pipefgmres\n")
   f.write("mpirun ../../maddm -N $N -problem ex2 -snes_type newtonls -ksp_rtol 1e-2 -snes_linesearch_order 2 -ksp_type pipefgmres\n")
   f.write("mpirun ../../maddm -N $N -problem ex3 -snes_type newtonls -ksp_rtol 1e-2 -snes_linesearch_order 2 -ksp_type pipefgmres\n")
   f.close()