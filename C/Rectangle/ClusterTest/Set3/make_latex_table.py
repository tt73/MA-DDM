import matplotlib.pyplot as plt
import numpy as np

# filename = 'out_test3'

# Reading files usually depends on 3 things:
# 1. Filename.
# 2. Number of trials. A trial consists of a run of each of the 3 problems.
# 3. Header size. This is the number of lines to skip for each trial.

print("Enter filename: ")
filename = input()
print("Enter number of trials: ")
Ntrials = int(input())
print("Enter header size: ")
Nheader = int(input())

f = open(filename)

N = Ntrials
variables = []
errors = np.zeros((3,N),dtype=float)
times = np.zeros((3,N),dtype=float)
iters = np.zeros((3,N),dtype=int)

for i in range(N):
   for j in range(Nheader):
      variables.append(f.readline()) # skip the header
   for j in range(3):
      f.readline() # Problem
      f.readline() # Params
      errors[j][i] = f.readline().split()[-1] # Error
      times[j][i] = f.readline().split()[-1]  # WTime
      iters[j][i] = f.readline().split()[-1]  # Iters

f.close()

print(variables)
print("Errors:")
print(errors)
print("Times:")
print(times)
print("Iters:")
print(iters)

print("LaTeX Table: ")
for i in range(3):
   print("P{:d} & ".format(i+1),end='')
   for j in range(N):
      print("{:.3f} & {:d}".format(times[i][j],iters[i][j]),end='')
      if (j<N-1):
         print(" & ",end='')
      else:
         print(" \\\\")
