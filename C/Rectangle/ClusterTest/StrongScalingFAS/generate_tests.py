import numpy as np
import subprocess

nps = np.array([2,4,6,8,9])

# loop to create job files
for i in range(len(nps)):
   jobfile = "job{:d}".format(nps[i])
   f = open(jobfile, "w")
   f.write("#!/bin/bash -l\n")
   f.write("#SBATCH -J j{:d}\n".format(nps[i]))
   f.write("#SBATCH -o slurmout{:d}\n".format(nps[i]))
   f.write("#SBATCH -p dms-cpu\n")
   f.write("#SBATCH -A tt73\n")
   f.write("#SBATCH -t 8:0:0\n")
   f.write("#SBATCH --nodes {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks-per-node 1\n")
   f.write("#SBATCH --exclusive\n")
   f.write("module load gnu8 mpich petsc/3.12.0\n")
   f.write("N=300\n")
   f.write("mpirun ../../maddm -N $N -problem ex1 -snes_type fas -op 0.1 >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../maddm -N $N -problem ex2 -snes_type fas -op 0.1 >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../maddm -N $N -problem ex3 -snes_type fas -op 0.1 >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../maddm -N $N -problem ex4 -snes_type fas -op 0.1 >> out{:d}\n".format(nps[i]))
   f.close()