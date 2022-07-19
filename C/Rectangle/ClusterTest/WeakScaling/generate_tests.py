import numpy as np
import subprocess

base_N = 100

nps = np.array([1,4,9])
Nxs = np.array([base_N, base_N*2, base_N*3])
Nys = np.array([base_N, base_N*2, base_N*3])

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
   f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -problem ex1 -sin >> out{:d}\n".format(Nxs[i],Nys[i],nps[i]))
   f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -problem ex2 -sin >> out{:d}\n".format(Nxs[i],Nys[i],nps[i]))
   f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -problem ex3 -sin >> out{:d}\n".format(Nxs[i],Nys[i],nps[i]))
   f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -problem ex4 -sin >> out{:d}\n".format(Nxs[i],Nys[i],nps[i]))
   f.close()