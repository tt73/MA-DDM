import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

nps = np.array([2,4,6,8,9])
n = len(nps)

# creates arrays for data
times = np.zeros((4,n),dtype=float)
errs = np.zeros((4,n),dtype=float)
iters = np.zeros((4,n),dtype=int)

# loop over all files

for i in range(n):
   f = open("out{}".format(nps[i]))
   # loop over each of the 3 problems
   # read in the error, runtime, and iterations
   for j in range(4):
      f.readline() # skip
      f.readline() # skip
      errs[j][i]  = f.readline().split()[-1]
      times[j][i] = f.readline().split()[-1]
      iters[j][i] = f.readline().split()[-1]
   f.close()


# base_times = np.array([22.089607, 7.997088, 36.583053]) # N = 200
base_times = np.array([97.899028, 29.104737, 176.157447, 38.895162]) # N = 300



# print(errs)
# print(times)
# print(iters)

# plot runtime vs Nd



plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",4)
fig = plt.figure(1, figsize=(6,6))

speedup1 = base_times[0]/times[0]
speedup2 = base_times[1]/times[1]
speedup3 = base_times[2]/times[2]
speedup4 = base_times[3]/times[3]

plt.plot(nps,speedup1,'-', label='P1',marker='o',mfc='w',color=colors[0],ms=8)
plt.plot(nps,speedup2,':', label='P2',marker='s',mfc='w',color=colors[1],ms=8)
plt.plot(nps,speedup3,'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
plt.plot(nps,speedup4,'.-',label='P4',marker='x',mfc='w',color=colors[3],ms=8)

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.legend(fontsize=14)
plt.xticks(nps)
plt.yticks(np.arange(1,10))
plt.xlabel('$N_d$')
plt.ylabel('Speedup')
plt.savefig('speedup.png',dpi=600,bbox_inches='tight')
plt.show()

