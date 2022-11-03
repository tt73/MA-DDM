import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
sys.path.append('../../')
import table_module as tm

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


# n = n_tols
# m = n_lims
# row_name = "Tols"
# col_name = "Limits"
# row_labs = tols
# col_labs = limits
# A = np.zeros((n,m),dtype=float)
# for i in range(n):
#    for j in range(m):
#       A[i,j] = iters[j,i,0]
# print(A)

# tm.print_table(n,m,row_name,col_name,row_labs,col_labs,A)

for i in range(6):
   for j in range(n_lims):
      print("{} {} ".format(iters[j,i,0],iters[j,i,1]),end="")
   print("")

for i in range(6):
   for j in range(n_lims):
      print("{} {} ".format(errs[j,i,0],errs[j,i,1]),end="")
   print("")

for i in range(6):
   for j in range(n_lims):
      print("{} {} ".format(times[j,i,0],times[j,i,1]),end="")
   print("")



plt.rcParams.update({'font.size' : 14})
# colors=sns.color_palette("rocket",4)
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.semilogx(tols,iters[0,:,0],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.semilogx(tols,iters[1,:,0],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.semilogx(tols,iters[2,:,0],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.semilogx(tols,iters[3,:,0],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

# # plt.yticks(np.arange(1,10))
# plt.ylim([0,1.1])
plt.xlabel('Tolerance')
plt.ylabel('Iterations')
plt.savefig('ex1_nd4_iters_h05.png',dpi=600,bbox_inches='tight')
# plt.show()
plt.clf()


plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.semilogx(tols,iters[0,:,1],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.semilogx(tols,iters[1,:,1],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.semilogx(tols,iters[2,:,1],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.semilogx(tols,iters[3,:,1],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

plt.xlabel('Tolerance')
plt.ylabel('Iterations')
plt.savefig('ex1_nd4_iters_h01.png',dpi=600,bbox_inches='tight')
plt.clf()



plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.loglog(tols,errs[0,:,0],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.loglog(tols,errs[1,:,0],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.loglog(tols,errs[2,:,0],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.loglog(tols,errs[3,:,0],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

plt.xlabel('Tolerance')
plt.ylabel('Error')
plt.savefig('ex1_nd4_errs_h05.png',dpi=600,bbox_inches='tight')
plt.clf()



plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.loglog(tols,errs[0,:,1],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.loglog(tols,errs[1,:,1],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.loglog(tols,errs[2,:,1],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.loglog(tols,errs[3,:,1],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

plt.xlabel('Tolerance')
plt.ylabel('Error')
plt.savefig('ex1_nd4_errs_h01.png',dpi=600,bbox_inches='tight')
plt.clf()



plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.loglog(tols,times[0,:,0],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.loglog(tols,times[1,:,0],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.loglog(tols,times[2,:,0],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.loglog(tols,times[3,:,0],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

plt.xlabel('Tolerance')
plt.ylabel('Runtime (s)')
plt.savefig('ex1_nd4_times_h05.png',dpi=600,bbox_inches='tight')
plt.clf()


plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("tab10",4)
fig = plt.figure(1, figsize=(6,6))

plt.loglog(tols,times[0,:,1],'-', color=colors[0],marker='o',ms=8,mfc='w',label='$L$=0.5')
plt.loglog(tols,times[1,:,1],':', color=colors[1],marker='s',ms=8,mfc='w',label='$L$=1.0')
plt.loglog(tols,times[2,:,1],'--',color=colors[2],marker='^',ms=8,mfc='w',label='$L$=1.5')
plt.loglog(tols,times[3,:,1],'.-',color=colors[3],marker='x',ms=8,mfc='w',label='$L$=2.0')

plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on()
plt.tick_params(which='minor', bottom=False, top=False)
plt.legend(fontsize=14)

plt.xlabel('Tolerance')
plt.ylabel('Runtime (s)')
plt.savefig('ex1_nd4_times_h01.png',dpi=600,bbox_inches='tight')
plt.clf()