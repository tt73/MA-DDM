# This is an example of running petsc programs.
import os
import subprocess as sp
import timeit 

# Function to run a bash command and then save the output to string `out`
# def run_petsc(command):
#    # p = sp.Popen(["/bin/bash", "-i","-c",command],stdout=sp.PIPE,universal_newlines=1,shell=1)
#    p = sp.check_output(command)
#    print(p)
#    return out


# The tri.c program simply solves a tridiagonal system of size m and then prints the error.
# Here is how it can be run
np = 2    # number of processors
m = 100000  # size of linear system
cmd = "mpiexec -n {} ./tri -tri_m {}".format(np,m)
print("The input:  "+cmd)
out = sp.check_output(cmd, shell=True, universal_newlines=True)
print("The output: "+out)

