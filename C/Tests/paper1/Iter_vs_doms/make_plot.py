import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

nps = np.array([4,9,16,25])
n = len(nps)
ops = np.array([0.1, 0.2, 0.3, 0.4])
m = len(ops)

# creates arrays for data
times = np.zeros((n,m,2),dtype=float)
errs = np.zeros((n,m,2),dtype=float)
iters = np.zeros((n,m,2),dtype=int)

# loop over all files
for i in range(n):
   f = open("out{}".format(nps[i]))
   for j in range(m):
      for k in range(2):
         f.readline() # skip
         f.readline() # skip
         f.readline() # skip
         f.readline() # skip
         errs [i][j][k]  = f.readline().split()[-1]
         times[i][j][k] = f.readline().split()[-1]
         iters[i][j][k] = f.readline().split()[-1]
   f.close()


print(errs)
print(times)
print(iters)

# plot runtime vs Nd



plt.rcParams.update({'font.size' : 20})
colors=sns.color_palette("rocket",2)
fig = plt.figure(1, figsize=(6,6))
ax = fig.add_axes([0,0,1,1])
ticklabels = ['2x2','3x3','4x4','5x5']

plt.plot(iters[:,1,0],'-', label='Ex1 p=20%',marker='s',mfc='w',color=colors[0],ms=8)
plt.plot(iters[:,1,1],'-', label='Ex2 p=20%',marker='o',mfc='w',color=colors[1],ms=8)
plt.plot(iters[:,3,0],'--',label='Ex1 p=40%',marker='s',mfc='w',color=colors[0],ms=8)
plt.plot(iters[:,3,1],'--',label='Ex2 p=40%',marker='o',mfc='w',color=colors[1],ms=8)

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=16)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.legend(fontsize=16)

plt.xticks(np.arange(4),ticklabels)

# plt.xticks(nps[1::])
# ax.set_xticklabels(labels)

# # plt.yticks(np.arange(0,10))
plt.xlabel('Domain Decomposition')
plt.ylabel('Iterations')
plt.savefig('test.png',dpi=100,bbox_inches='tight')
# plt.show()

