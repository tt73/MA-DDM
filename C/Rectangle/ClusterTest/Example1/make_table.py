import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


limits = np.array([0.5, 1.0, 1.5, 2.0])
n_lims = len(limits)

tols   = np.array([1e-1, 1e-4, 1e-6])
n_tols = len(tols)

sizes  = np.arrary([100, 200, 300, 400])
n_sizes = len(sizes)


# creates arrays for data
times = np.zeros((n_lims,n_tols,n_sizes),dtype=float)
errs = np.zeros((n_lims,n_tols,n_sizes),dtype=float)
iters = np.zeros((n_lims,n_tols,n_sizes),dtype=int)

# choose file
file = "Nd1.out"
f = open(file)
for i in range(n_lims):
   for ii in range(n_tols):
      for iii in range(n_sizes):
         f.readline() # skip
         f.readline() # skip
         errs[i][ii][iii]  = f.readline().split()[-1]
         times[i][ii][iii] = f.readline().split()[-1]
         iters[i][ii][iii] = f.readline().split()[-1]
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
plt.ylim([0,1.1])
plt.xlabel('$N_d$')
plt.ylabel('Efficiency')
plt.savefig('efficiency_htn.png',dpi=600,bbox_inches='tight')
plt.show()