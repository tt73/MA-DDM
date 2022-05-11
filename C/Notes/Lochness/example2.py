# This is an example of running petsc programs.
import subprocess as sp
import time

# Function to run a bash command and then save the output to string `out`
def run_tri(np,m):
   cmd = "/usr/bin/mpiexec -n {} ./tri -tri_m {}".format(np,m)
   # cmd = "/usr/bin/mpiexec -n {} ./tri -tri_m {}".format(np,m)
   p = sp.run([cmd],capture_output=1,universal_newlines=1,shell=1)
   return p


# The tri.c program simply solves a tridiagonal system of size m and then prints the error.
# Here is how it can be run:
np = 2    # number of processors
m = 100000  # size of linear system
res1 = run_tri(np,m)
print(res1.stdout)


# You can time the code with the Time module:
m = 1000000
runtimes = []
for i in range(8):
   tic = time.time()
   run_tri(i+1,m)
   toc = time.time()
   print("n = {:3d} | time = {:10.4f}".format(i+1,toc-tic))
