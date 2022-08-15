import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",4)

f = open("out1")
N = 4
M = 4
times = np.zeros((N,M,4),dtype=float)
errs = np.zeros((N,M,4),dtype=float)
iters = np.zeros((N,M,4),dtype=int)

for i in range(N):
   for j in range(M):
      print(f.readline()) # tolerance
      for k in range(4):
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

print("problem 4")
print(times[:,:,3])
print(iters[:,:,3])

_x = np.arange(N)
_y = np.arange(M)
_xx, _yy = np.meshgrid(_x, _y)
x, y = _xx.ravel(), _yy.ravel()
z0 = np.zeros_like(x)
dx = 1
dy = 1


fig = plt.figure(figsize=(6, 6))
ax1=fig.add_subplot(111, projection='3d')
ax1.set_xlabel('Newton tol', labelpad=10)
ax1.set_ylabel('Krylov tol', labelpad=10)
ax1.set_zlabel('Runtime (sec)')
ax1.view_init(20,51)

ax1.bar3d(x,y,z0,dx,dy,times[:,:,0].ravel(),color=colors[0],alpha=0.8)
names = ['1e-4','1e-3','1e-2','1e-1']
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN1.png',dpi=600,bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(6, 6))
ax2=fig.add_subplot(111, projection='3d')
ax2.set_xlabel('Newton tol', labelpad=10)
ax2.set_ylabel('Krylov tol', labelpad=10)
ax2.set_zlabel('Runtime (sec)')
ax2.view_init(20,51)
ax2.bar3d(x,y,z0,dx,dy,times[:,:,1].ravel(),color=colors[1],alpha=0.8)
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN2.png',dpi=600,bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(6, 6))
ax3=fig.add_subplot(111, projection='3d')
ax3.set_xlabel('Newton tol', labelpad=10)
ax3.set_ylabel('Krylov tol', labelpad=10)
ax3.set_zlabel('Runtime (sec)')
ax3.view_init(20,51)
ax3.bar3d(x,y,z0,dx,dy,times[:,:,2].ravel(),color=colors[2],alpha=0.8)
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN3.png',dpi=600,bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(6, 6))
ax3b=fig.add_subplot(111, projection='3d')
ax3b.set_xlabel('Newton tol', labelpad=10)
ax3b.set_ylabel('Krylov tol', labelpad=10)
ax3b.set_zlabel('Iters')
ax3b.view_init(20,51)
ax3b.bar3d(x,y,z0,dx,dy,iters[:,:,2].ravel(),color=colors[2],alpha=0.8)
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN3_iter.png',dpi=600,bbox_inches='tight')
plt.show()

fig = plt.figure(figsize=(6, 6))
ax3b=fig.add_subplot(111, projection='3d')
ax3b.set_xlabel('Newton tol', labelpad=10)
ax3b.set_ylabel('Krylov tol', labelpad=10)
ax3b.set_zlabel('h-norm Error')
ax3b.view_init(20,51)
ax3b.bar3d(x,y,z0,dx,dy,errs[:,:,2].ravel(),color=colors[2],alpha=0.8)
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN3_err.png',dpi=600,bbox_inches='tight')
plt.show()


fig = plt.figure(figsize=(6, 6))
ax4=fig.add_subplot(111, projection='3d')
ax4.set_xlabel('Newton tol', labelpad=10)
ax4.set_ylabel('Krylov tol', labelpad=10)
ax4.set_zlabel('Runtime (sec)')
ax4.view_init(20,51)
ax4.bar3d(x,y,z0,dx,dy,times[:,:,3].ravel(),color=colors[3],alpha=0.8)
plt.xticks(_x,names)
plt.yticks(_y,names)
plt.savefig('HTN4.png',dpi=600,bbox_inches='tight')
plt.show()