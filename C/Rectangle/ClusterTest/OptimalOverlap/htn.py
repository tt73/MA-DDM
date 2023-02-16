import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

sys.path.append('../../')
import table_module as tm

# choose file
file = "sweephtn.out"

dim1 = 4 # newton
dim2 = 4 # krylov
dim3 = 3 # overlap

# creates arrays for data
times = np.zeros((dim1,dim2,dim3),dtype=float)
errs = np.zeros((dim1,dim2,dim3),dtype=float)
iters = np.zeros((dim1,dim2,dim3),dtype=int)

f = open(file)
for i in range(dim1):
   for ii in range(dim2):
      for iii in range(dim3):
         f.readline()
         f.readline()
         f.readline()
         f.readline()
         errs[i][ii][iii]  = f.readline().split()[-1]
         times[i][ii][iii] = f.readline().split()[-1]
         iters[i][ii][iii] = f.readline().split()[-1]
f.close()


print()
print(times)
print()
print(iters)
print()
print(errs)

print(np.amin(times))
print(np.argmin(times))

# for i in range(dim1):
#    for j in range(dim2):

#       print('{:8.4e} '.format(times[i][j]),end='')
#    print('')


# # h = 0.05
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,0],errs[i,-1,0],times[i,-1,0]))

# # h = 0.01
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,1],errs[i,-1,1],times[i,-1,1]))
