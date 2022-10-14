import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# choose file
file = "Nd4.out"
# file = "Nd4.out"

if (file=="Nd9.out"):
   limits = np.array([0.5, 1.0, 1.5])
else:
   limits = np.array([0.5, 1.0, 1.5, 2.0])
n_lims = len(limits)

tols   = np.array([1e-1,1e-2,1e-3,1e-4,1e-5,1e-6])
n_tols = len(tols)
print(n_tols)

if (file=="Nd1.out"):
   sizes = np.array([0.05, 0.01, 0.005])
else:
   sizes  = np.array([0.05, 0.01])
n_sizes = len(sizes)

# creates arrays for data
times = np.zeros((n_lims,n_tols,n_sizes),dtype=float)
errs = np.zeros((n_lims,n_tols,n_sizes),dtype=float)
iters = np.zeros((n_lims,n_tols,n_sizes),dtype=int)

f = open(file)
for i in range(n_lims):
   for ii in range(n_tols): # 6 tols
      for iii in range(n_sizes): # 2 hs
         f.readline() # skip
         f.readline() # skip
         if (file=='Nd1.out'):
            f.readline() # skip extra line for -snes_converged_reason
         errs[i][ii][iii]  = f.readline().split()[-1]
         times[i][ii][iii] = f.readline().split()[-1]
         iters[i][ii][iii] = f.readline().split()[-1]
f.close()

# print(errs)
# print(times)


print("cols: h = 0.05, 0.01")
print("rows:\ntol=\n 1e-1,\n 1e-2,\n,1e-3,\n1e-4,\n1e-5,\n 1e-6\n")
for i in range(n_lims):
   print("\nDomain [-{:.2f}, {:.2f}]^2\n".format(limits[i],limits[i]))
   print("Iters:")
   print(iters[i,:,:])
   print("Error:")
   print(errs[i,:,:])
   print("Time:")
   print(times[i,:,:])


for j in range(n_lims):
   print("{} {} ".format("h=0.05","h=0.01"),end="")
print("")


# for i in range(6):
#    for j in range(n_lims):
#       print("{} {} ".format(iters[j,i,0],iters[j,i,1]),end="")
#    print("")

# for i in range(6):
#    for j in range(n_lims):
#       print("{} {} ".format(errs[j,i,0],errs[j,i,1]),end="")
#    print("")

for i in range(6):
   for j in range(n_lims):
      print("{} {} ".format(times[j,i,0],times[j,i,1]),end="")
   print("")

# plot runtime vs Nd

# plt.rcParams.update({'font.size' : 14})
# colors=sns.color_palette("rocket",4)
# fig = plt.figure(1, figsize=(6,6))

# plt.plot(nps,times[0,0]/times[0],'-', label='P1',marker='o',mfc='w',color=colors[0],ms=8)
# plt.plot(nps,times[1,0]/times[1],':', label='P2',marker='s',mfc='w',color=colors[1],ms=8)
# plt.plot(nps,times[2,0]/times[2],'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
# plt.plot(nps,times[3,0]/times[3],'.-',label='P4',marker='x',mfc='w',color=colors[3],ms=8)

# plt.tick_params(labelsize=14)
# plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
# plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
# plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
# plt.minorticks_on()
# plt.tick_params(which='minor', bottom=False, top=False)
# plt.legend(fontsize=14)
# plt.xticks(nps)


# # plt.yticks(np.arange(1,10))
# plt.ylim([0,1.1])
# plt.xlabel('$N_d$')
# plt.ylabel('Efficiency')
# plt.savefig('efficiency_htn.png',dpi=600,bbox_inches='tight')
# plt.show()