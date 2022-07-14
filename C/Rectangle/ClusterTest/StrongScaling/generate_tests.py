import numpy as np
import subprocess

nps = np.array([2,4,6,8,9])

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
   # f.write("#SBATCH --exclusive\n")
   f.write("module load gnu8 mpich petsc\n")
   f.write("N=200\n")
   f.write("mpirun ../../maddm -N $N -problem ex1 -sin >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../maddm -N $N -problem ex2 -sin >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../maddm -N $N -problem ex3 -sin >> out{:d}\n".format(nps[i]))
   f.close()