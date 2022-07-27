import numpy as np
import subprocess

base_N = 100
op = 0.3 # overlap percentage

nps = np.array([1,4,9])
Nxs = np.array([base_N, base_N*2/(1+op), base_N*3/(1+op)],dtype=int)
Nys = np.array([base_N, base_N*2/(1+op), base_N*3/(1+op)],dtype=int)

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
   f.write("#SBATCH --mem=0G\n")
   f.write("#SBATCH --nodes {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks {:d}\n".format(nps[i]))
   f.write("module load gnu8 mpich petsc\n")
   if (i==0):
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex1 >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex2 >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex3 >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex4 >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))

   else:
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex1 -htn >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex2 -htn >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex3 -htn >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
      f.write("mpirun ../../maddm -Nx {:d} -Ny {:d} -op {:f} -problem ex4 -htn >> out{:d}\n".format(Nxs[i],Nys[i],op,nps[i]))
   f.close()