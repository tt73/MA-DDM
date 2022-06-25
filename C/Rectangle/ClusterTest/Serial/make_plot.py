import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

colors=sns.color_palette("rocket",3)

f = open("out_test1b")
N = 9
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

fig = plt.figure(1, figsize=(6,6))
plt.plot(Ns,times[0],'-', label='P1',marker='o',mfc='w',color=colors[0])
plt.plot(Ns,times[1],':', label='P2',marker='s',mfc='w',color=colors[1])
plt.plot(Ns,times[2],'--',label='P3',marker='^',mfc='w',color=colors[2])
plt.xlabel('$N$')
plt.ylabel('Time (sec)')
plt.savefig('test.png')
plt.legend()
plt.show()

