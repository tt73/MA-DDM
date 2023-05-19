import numpy as np
import subprocess


# loop to create parallel job files
nps = np.array([2,4,9,16,25,36,49])
for i in range(len(nps)):
   jobfile = "job{:d}".format(nps[i])
   f = open(jobfile, "w")
   f.write("#!/bin/bash -l\n")
   f.write("#SBATCH -J j{:d}\n".format(nps[i]))
   f.write("#SBATCH -o slurmout{:d}\n".format(nps[i]))
   f.write("#SBATCH -p lowpriority\n")
   f.write("#SBATCH --constraint=avx512\n")
   f.write("#SBATCH -A tt73\n")
   f.write("#SBATCH -t 8:0:0\n")
   f.write("#SBATCH --mem=0G\n")
   f.write("#SBATCH --nodes {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks {:d}\n".format(nps[i]))
   f.write("#SBATCH --ntasks-per-node 1\n")
   f.write("module load gnu8 mpich petsc/3.12.0\n")
   f.write("N=401\n")
   f.write("for op in {0.1,0.2,0.3,0.4}\n")
   f.write("do\n")
   f.write("   mpirun ../../maddm -N $N -problem ex1 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -init_type coarse -coarseness 4 -snes_converged_reason -op $op >> out{:d}\n".format(nps[i]))
   f.write("   mpirun ../../maddm -N $N -problem ex5 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -init_type coarse -coarseness 4 -snes_converged_reason -op $op >> out{:d}\n".format(nps[i]))
   f.write("done")
   f.close()


