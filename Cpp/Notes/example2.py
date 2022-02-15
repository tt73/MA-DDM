# This is an example of running petsc programs.

import subprocess as sp


# Function to run a bash command and then save the output to string `out`
def run_petsc(command):
   p = sp.Popen(["/bin/bash", "-i","-c",command],stdout=sp.PIPE,universal_newlines=1,shell=1)
   out,err = p.communicate()
   return out


# The tri.c program simply solves a tridiagonal system of size m and then prints the error.
# Here is how it can be run
np = 2    # number of processors
m = 1000  # size of linear system
cmd = "petsc -n {} ./tri -tri_m {}".format(np,m)

out = run_petsc(cmd)
print(out)


# You can time the code with the command time
cmd = "time "+cmd
print(cmd)
out = run_petsc(cmd)
print(out)