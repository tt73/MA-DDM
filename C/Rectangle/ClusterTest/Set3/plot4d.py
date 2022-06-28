import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",3)

f = open("out4d")
N = 4
M = 1
times = np.zeros((N,M,3),dtype=float)
errs = np.zeros((N,M,3),dtype=float)
iters = np.zeros((N,M,3),dtype=int)

for i in range(N):
   for j in range(M):
      print(f.readline()) # tolerance
      for k in range(3):
         f.readline() # Problem
         f.readline() # Params
         # print(f.readline().split())
         errs[i][j][k] = f.readline().split()[-1]  # Error
         times[i][j][k] = f.readline().split()[-1] # WTime
         iters[i][j][k] = f.readline().split()[-1] # Iters

print(errs)

print("problem 1")
print(times[:,:,0])
print(iters[:,:,0])

print("problem 2")
print(times[:,:,1])
print(iters[:,:,1])

print("problem 3")
print(times[:,:,2])
print(iters[:,:,2])

fig = plt.figure(figsize=(6, 6))
# ax1=fig.add_subplot(111, projection='3d')
ax1=fig.add_subplot(111)
ax1.set_xlabel('rtol', labelpad=10)
ax1.set_ylabel('Time (sec)', labelpad=10)

x = np.arange(4)

ax1.plot(x,times[:,:,0].ravel(),'-',label='P1',marker='o',mfc='w',color=colors[0],ms=8)
ax1.plot(x,times[:,:,1].ravel(),':',label='P2',marker='s',mfc='w',color=colors[1],ms=8)
ax1.plot(x,times[:,:,2].ravel(),'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
names = ['1e-4','1e-3','1e-2','1e-1']
plt.xticks(x,names)
# plt.yticks(_y,names)
# plt.savefig('rtols.png',dpi=300,bbox_inches='tight')
plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.legend(fontsize=14)
plt.savefig('fig3.png',dpi=300,bbox_inches='tight')
plt.show()




