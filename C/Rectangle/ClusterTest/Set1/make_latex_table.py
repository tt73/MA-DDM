import matplotlib.pyplot as plt
import numpy as np

filename = 'out_test3'

f = open(filename)
N = 6
op = np.zeros(N,dtype=float)
times = np.zeros((3,N),dtype=float)
iters = np.zeros((3,N),dtype=int)
for i in range(N):
   f.readline() # N
   for j in range(3):
      f.readline() # Problem
      f.readline() # Params
      f.readline() # Error
      times[j][i] = f.readline().split()[1] # WTime
      iters[j][i] = f.readline().split()[1] # Iters
   op[i] = i*0.05

print(times)
print(iters)

for i in range(3):
   print("P{:d} & ".format(i),end='')
   for j in range(N):
      print("{:.3f} & {:d}".format(times[i][j],iters[i][j]),end='')
      if (j<N-1):
         print(" & ",end='')
      else:
         print(" \\\\")