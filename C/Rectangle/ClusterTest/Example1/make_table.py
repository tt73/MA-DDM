import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

sys.path.append('../../')
import table_module as tm

# choose file
Nd = 4
loc = 2
glo = 9
file = "Nd{}_loc{}_glo{}.out".format(Nd,loc,glo)

# file = "Nd9.out"

if (file=="Nd6.out" or file=="Nd8.out" or file=="Nd9.out"):
   limits = np.array([1.5, 2.0])
else:
   limits = np.array([0.5, 1.0, 1.5, 2.0])
n_lims = len(limits)

if (file=="Nd6.out" or file=="Nd8.out" or file=="Nd9.out"):
# if (file=="Nd4.out" or file=="Nd6.out" or file=="Nd8.out" or file=="Nd9.out"):
   tols = np.array([1e-2])
elif(file=="Nd4.out"):
   tols = np.array([1e-2])
   # tols = np.array([1e-1,1e-2])
else:
   tols = np.array([1e-1,1e-2,1e-3,1e-4,1e-5,1e-6])
n_tols = len(tols)

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
         f.readline() # skip extra line for -snes_converged_reason
         f.readline() # skip
         f.readline() # skip
         errs[i][ii][iii]  = f.readline().split()[-1]
         times[i][ii][iii] = f.readline().split()[-1]
         iters[i][ii][iii] = f.readline().split()[-1]
         print(iters[i][ii][iii])
f.close()



# print(times)

print(errs)


# h = 0.05
print("limit, iter, error, runtime")
for i in range(n_lims):
   print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,0],errs[i,-1,0],times[i,-1,0]))

# h = 0.01
print("limit, iter, error, runtime")
for i in range(n_lims):
   print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,1],errs[i,-1,1],times[i,-1,1]))


# # h = 0.05
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,0,0],errs[i,0,0],times[i,0,0]))

# # h = 0.01
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,0,1],errs[i,0,1],times[i,0,1]))

