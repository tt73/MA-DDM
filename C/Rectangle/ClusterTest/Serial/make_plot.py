import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

f = open("out1")
N = 9
Ns = np.zeros(N,dtype=int)
times = np.zeros((4,N),dtype=float)
for i in range(N):
   f.readline() # N
   for j in range(4):
      f.readline() # Problem
      f.readline() # Params
      f.readline() # Error
      times[j][i] = f.readline().split()[1] # WTime
      f.readline() # Iters
   Ns[i] = 100 + i*50

plt.rcParams.update({'font.size' : 14})
colors = sns.color_palette("rocket",4)
lines = ['-',':','--','.-']
marks = ['o','s','^','x']
fig = plt.figure(1, figsize=(6,6))
plt.plot(Ns,times[0],lines[0],label='P1',marker=marks[0],mfc='w',color=colors[0],ms=8)
plt.plot(Ns,times[1],lines[1],label='P2',marker=marks[1],mfc='w',color=colors[1],ms=8)
plt.plot(Ns,times[2],lines[2],label='P3',marker=marks[2],mfc='w',color=colors[2],ms=8)
plt.plot(Ns,times[3],lines[3],label='P4',marker=marks[3],mfc='w',color=colors[3],ms=8)
plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.xlabel('$N$')
plt.ylabel('Runtime (sec)')
plt.yticks(np.arange(0,1201,120))
plt.legend(fontsize=14)
plt.savefig('serial.png',dpi=600,bbox_inches='tight')
plt.show()

