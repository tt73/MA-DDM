import matplotlib.pyplot as plt
import numpy as np
import sys

f = open("out_test1")
N = 7
Ns = np.zeros(N,dtype=int)
times = np.zeros((3,N),dtype=float)
for i in range(N):
   f.readline() # N
   for j in range(3):
      f.readline() # Problem
      f.readline() # Params
      f.readline() # Error
      times[j][i] = f.readline().split()[1] # WTime
      f.readline() # Iters
   Ns[i] = 100 + i*50

fig = plt.figure(1, figsize=(5,5))
plt.plot(Ns,times[0],'-', label='P1',marker='o',mfc='w')
plt.plot(Ns,times[1],':', label='P2',marker='s',mfc='w')
plt.plot(Ns,times[2],'--',label='P3',marker='^',mfc='w')
plt.xlabel('$N$')
plt.ylabel('Time (sec)')
plt.savefig('test.png')
plt.legend()
plt.show()

