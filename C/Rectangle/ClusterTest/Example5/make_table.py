import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# data settings
dim1 = 4
dim2 = 2
dim3 = 3
Nheader = 1
fname  = 'out_L'
var_labels = ["HTN","SIN","NKS"]

f = open(fname)
errors = np.zeros((dim1,dim2,dim3),dtype=float)
times  = np.zeros((dim1,dim2,dim3),dtype=float)
iters  = np.zeros((dim1,dim2,dim3),dtype=int)

# Loop over all trials
#    Skip the header lines
#    Loop over the 4 examples
for i in range(dim1):
   for j in range(dim2):
      for k in range(dim3):
         f.readline() # converged reason
         f.readline() # Problem
         f.readline() # Params
         errors[i][j][k] = f.readline().split()[-1] # Error
         times[i][j][k] = f.readline().split()[-1]  # WTime
         iters[i][j][k] = f.readline().split()[-1]  # Iters
f.close()

## print out tables

print("runtime")
rowstr = ""
for i in range(dim1):
   rowstr += "{:^28f}".format(0.5*(i+1))
print(rowstr)
for k in range(dim3):
   rowstr = ""
   for i in range(dim1):
      for j in range(dim2):
         rowstr += "{:12.3f} &".format(times[i,j,k])
   print(rowstr)

print("error")
rowstr = ""
for i in range(dim1):
   rowstr += "{:^28f}".format(0.5*(i+1))
print(rowstr)
for k in range(dim3):
   rowstr = ""
   for i in range(dim1):
      for j in range(dim2):
         rowstr += "{:12.3E} &".format(errors[i,j,k])
   print(rowstr)
# rowstr = ""
# for i in range(Ntrials):
#    rowstr += "{:>12s} &".format(variables[i])
# print(rowstr)

# for j in range(Nvariables):
#    rowstr = ""
#    for i in range(Ntrials):
#       rowstr += "{:12.3E} &".format(errors[j,i])
#    print(rowstr)

print("iterations")
rowstr = ""
for i in range(dim1):
   rowstr += "{:^28f}".format(0.5*(i+1))
print(rowstr)
for k in range(dim3):
   rowstr = ""
   for i in range(dim1):
      for j in range(dim2):
         rowstr += "{:12d} &".format(iters[i,j,k])
   print(rowstr)
# rowstr = ""
# for i in range(Ntrials):
#    rowstr += "{:>12s} &".format(variables[i])
# print(rowstr)
# for j in range(Nvariables):
#    rowstr = ""
#    for i in range(Ntrials):
#       rowstr += "{:12d} &".format(iters[j,i])
#    print(rowstr)
