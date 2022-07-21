import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

nps = np.array([1,4,9])
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


print(errs)
print(times)
print(iters)

# plot runtime vs Nd

plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",4)
fig = plt.figure(1, figsize=(6,6))

plt.plot(nps,times[0,0]/times[0],'-', label='P1',marker='o',mfc='w',color=colors[0],ms=8)
plt.plot(nps,times[1,0]/times[1],':', label='P2',marker='s',mfc='w',color=colors[1],ms=8)
plt.plot(nps,times[2,0]/times[2],'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
plt.plot(nps,times[3,0]/times[3],'.-',label='P4',marker='x',mfc='w',color=colors[3],ms=8)

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)
plt.xticks(nps)


# plt.yticks(np.arange(1,10))
plt.xlabel('$N_d$')
plt.ylabel('Efficiency')
plt.savefig('Efficiency.png',dpi=600,bbox_inches='tight')
plt.show()