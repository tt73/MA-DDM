import numpy as np
import subprocess

nps = np.array([2,4,9])

# loop to create job files
for i in range(len(nps)):
   jobfile = "job{:d}".format(nps[i])
   f = open(jobfile, "w")
   f.write("#!/bin/bash -l\n")
   f.write("#SBATCH -J j{:d}\n".format(nps[i]))
   f.write("#SBATCH -p public\n")
   f.write("#SBATCH -o slurmout{:d}\n".format(nps[i]))
   f.write("#SBATCH -p dms-cpu\n")
   f.write("#SBATCH --mail-user tt73@njit.edu\n")
   f.write("#SBATCH -A tt73\n")
   f.write("#SBATCH -t 8:0:0\n")
   f.write("#SBATCH --mem=16G\n")
   f.write("#SBATCH --nodes {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks-per-node 1\n")
   f.write("module load gnu8 mpich petsc/3.12.0\n")

   f.write("N=450\n")
   f.write("rm -f out{:d}\n".format(nps[i]))
   f.write("mpiexec -np {:d} ../../maddm -N $N -problem ex5 -htn -op 0.10 -sub_snes_rtol 1e-2 >> out{:d}\n".format(nps[i],nps[i]))
   # f.write("mpiexec -np {:d} ../../maddm -N $N -problem ex5 -sin -op 0.10 >> out{:d}\n".format(nps[i],nps[i]))
   f.write("mpiexec -np {:d} ../../maddm -N $N -problem ex5 -nks -op 0 -ksp_rtol 1e-2 >> out{:d}\n".format(nps[i],nps[i]))
   f.close()