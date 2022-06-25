import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

colors=sns.color_palette("rocket",3)

filename = 'out_test1'

f = open(filename)
N = 5
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
      print("{:.3f} & {:d} & ".format(times[i][j],iters[i][j]),end='')
   print('')




# fig = plt.figure(1, figsize=(6,6))
# plt.plot(Ns,times[0],'-', label='P1',marker='o',mfc='w',color=colors[0])
# plt.plot(Ns,times[1],':', label='P2',marker='s',mfc='w',color=colors[1])
# plt.plot(Ns,times[2],'--',label='P3',marker='^',mfc='w',color=colors[2])
# plt.tick_params(direction='in',right=True,top=True)
# plt.tick_params(labelsize=14)
# plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
# plt.xlabel('$N$')
# plt.ylabel('Time (sec)')
# plt.legend(fontsize=14)
# plt.savefig('serial.png',dpi=300,bbox_inches='tight')
# plt.show()

