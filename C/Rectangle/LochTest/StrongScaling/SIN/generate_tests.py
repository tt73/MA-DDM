import numpy as np
import subprocess



# serial code
jobfile = "job{:d}".format(1)
f = open(jobfile, "w")
f.write("#!/bin/bash -l\n")
f.write("#SBATCH -J j{:d}\n".format(1))
f.write("#SBATCH -o slurmout{:d}\n".format(1))
f.write("#SBATCH -p lowpriority\n")
f.write("#SBATCH --constraint=avx512\n")
f.write("#SBATCH -A tt73\n")
f.write("#SBATCH -t 8:0:0\n")
f.write("#SBATCH --mem=0G\n")
f.write("#SBATCH --nodes {:d}\n".format(1))
f.write("#SBATCH --ntasks {:d}\n".format(1))
f.write("#SBATCH --ntasks-per-node 1\n")
f.write("module load gnu8 mpich petsc/3.12.0\n")
f.write("N=300\n")
f.write("../../../maddm -N $N -problem ex1 >> out{:d}\n".format(1))
f.write("../../../maddm -N $N -problem ex2 >> out{:d}\n".format(1))
f.write("../../../maddm -N $N -problem ex3 >> out{:d}\n".format(1))
f.write("../../../maddm -N $N -problem ex4 >> out{:d}\n".format(1))
f.close()


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
   f.write("N=300\n")
   f.write("op=0.10\n")
   f.write("mpirun ../../../maddm -N $N -problem ex1 -sin -op $op >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../../maddm -N $N -problem ex2 -sin -op $op >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../../maddm -N $N -problem ex3 -sin -op $op >> out{:d}\n".format(nps[i]))
   f.write("mpirun ../../../maddm -N $N -problem ex4 -sin -op $op >> out{:d}\n".format(nps[i]))
   f.close()