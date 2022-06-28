import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",3)

f = open("out4a")
N = 4
M = 4
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
ax1=fig.add_subplot(111, projection='3d')
ax1.set_xlabel('SNES rtol', labelpad=10)
ax1.set_ylabel('KSP rtol', labelpad=10)
ax1.set_zlabel('Runtime')

_x = np.arange(N)
_y = np.arange(M)
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()
z0 = np.zeros_like(x)
dx = 1
dy = 1

# ax1.bar3d(x,y,z0,dx,dy,times[:,:,0].ravel(),color=colors[0])
# ax1.bar3d(x+0.2,y,z0,dx,dy,times[:,:,1].ravel(),color=colors[1])
ax1.bar3d(x+0.4,y,z0,dx,dy,times[:,:,2].ravel(),color=colors[2],alpha=0.8)
names = ['1e-4','1e-3','1e-2','1e-1']
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('rtols.png',dpi=300,bbox_inches='tight')
plt.show()


# fig = plt.figure(1, figsize=(6,6))
# plt.plot(Ns,times[0],'-', label='P1',marker='o',mfc='w',color=colors[0],ms=8)
# plt.plot(Ns,times[1],':', label='P2',marker='s',mfc='w',color=colors[1],ms=8)
# plt.plot(Ns,times[2],'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
# plt.tick_params(direction='in',right=True,top=True)
# plt.tick_params(labelsize=14)
# plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
# plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
# plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
# plt.xlabel('$N$')
# plt.ylabel('Time (sec)')
# plt.yticks(np.arange(0,1200,120))
# plt.legend(fontsize=14)
# plt.show()


