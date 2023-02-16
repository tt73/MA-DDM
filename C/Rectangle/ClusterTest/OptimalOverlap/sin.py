import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

sys.path.append('../../')
import table_module as tm

# choose file
file = "sweepa.out"

dim1 = 2 #
dim2 = 8 #

# creates arrays for data
times = np.zeros((dim1,dim2),dtype=float)
errs = np.zeros((dim1,dim2),dtype=float)
iters = np.zeros((dim1,dim2),dtype=int)

f = open(file)
for i in range(dim1):
   for ii in range(dim2):
      f.readline()
      f.readline()
      f.readline()
      f.readline()
      errs[i][ii]  = f.readline().split()[-1]
      times[i][ii] = f.readline().split()[-1]
      iters[i][ii] = f.readline().split()[-1]
f.close()


print()
print(times)
print()
print(iters)
print()
print(errs)

for i in range(dim1):
   for j in range(dim2):
      print('{:8.4e} '.format(times[i][j]),end='')
   print('')


# # h = 0.05
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,0],errs[i,-1,0],times[i,-1,0]))

# # h = 0.01
# print("limit, iter, error, runtime")
# for i in range(n_lims):
#    print("{} & {:3d} & {:8.6f} & {:8.4f} \\\\".format(limits[i],iters[i,-1,1],errs[i,-1,1],times[i,-1,1]))
