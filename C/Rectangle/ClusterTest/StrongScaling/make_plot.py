import numpy as np

# loop over all files
nps = np.array([2,4,6,8,9])
n = len(nps)
for i in range(n):
   f = open("out{}".format(nps[i]))
   # loop over each of the 3 problems
   # read in the error, runtime, and iterations
   for i in range(3):
      err = f.readline()
      time = f.readline()
      iter = f.readline()
      print(err)
      print(time)
      print(iter)
   f.close()
   print("done with file ")


base_times = np.array([12.2 10.1 123.12])


# plot runtime vs Nd

# import matplotlib.pyplot as plt
# import numpy as np
# import seaborn as sns

# plt.rcParams.update({'font.size' : 14})
# colors=sns.color_palette("rocket",3)

# f = open("out_test1b")
# N = 9
# Ns = np.zeros(N,dtype=int)
# times = np.zeros((3,N),dtype=float)
# for i in range(N):
#    f.readline() # N
#    for j in range(3):
#       f.readline() # Problem
#       f.readline() # Params
#       f.readline() # Error
#       times[j][i] = f.readline().split()[1] # WTime
#       f.readline() # Iters
#    Ns[i] = 100 + i*50

# fig = plt.figure(1, figsize=(6,6))
# plt.plot(Ns,times[0],'-', label='P1',marker='o',mfc='w',color=colors[0],ms=8)
# plt.plot(Ns,times[1],':', label='P2',marker='s',mfc='w',color=colors[1],ms=8)
# plt.plot(Ns,times[2],'--',label='P3',marker='^',mfc='w',color=colors[2],ms=8)
# plt.tick_params(direction='in',right=True,top=True)
# plt.tick_params(labelsize=14)
# plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
# plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
# plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
# plt.xlabel('$N$')
# plt.ylabel('Time (sec)')
# plt.yticks(np.arange(0,1200,120))
# plt.legend(fontsize=14)
# plt.savefig('serial.png',dpi=300,bbox_inches='tight')
# plt.show()

