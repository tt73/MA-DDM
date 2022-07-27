import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

Ntrials = 4
Nheader = 1

# read SIN
f = open('out2a')
N = Ntrials
variables = []
errors = np.zeros((4,N),dtype=float)
times = np.zeros((4,N),dtype=float)
iters = np.zeros((4,N),dtype=int)

for i in range(N):
   for j in range(Nheader):
      variables.append(f.readline().split()[-1]) # skip the header
   for j in range(4):
      f.readline() # Problem
      f.readline() # Params
      errors[j][i] = f.readline().split()[-1] # Error
      times[j][i] = f.readline().split()[-1]  # WTime
      iters[j][i] = f.readline().split()[-1]  # Iters
f.close()



plt.rcParams.update({'font.size' : 14})
colors = sns.color_palette("rocket",4)
lines = ['-',':','--','.-']
marks = ['o','s','^','x']


## Runtime for SIN
fig = plt.figure(1, figsize=(6,6))
plt.plot(variables,times[0],lines[0],label='P1',marker=marks[0],mfc='w',color=colors[0],ms=8)
plt.plot(variables,times[1],lines[1],label='P2',marker=marks[1],mfc='w',color=colors[1],ms=8)
plt.plot(variables,times[2],lines[2],label='P3',marker=marks[2],mfc='w',color=colors[2],ms=8)
plt.plot(variables,times[3],lines[3],label='P4',marker=marks[3],mfc='w',color=colors[3],ms=8)

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)

plt.xlabel('N')
plt.ylabel('Runtime (sec)')
# plt.yticks(np.arange(0,111,10))
plt.ylim([0,225])
plt.xticks(variables)

plt.legend(fontsize=14)
plt.savefig('time_sin2.png',dpi=600,bbox_inches='tight')
plt.show()


## Iteartions for SIN
fig = plt.figure(2, figsize=(6,6))
ax = fig.add_axes([0,0,1,1])
X = np.arange(N)
ax.bar(X + 0.0, iters[0], color = colors[0], width = 0.2, label='P1')
ax.bar(X + 0.2, iters[1], color = colors[1], width = 0.2, label='P2')
ax.bar(X + 0.4, iters[2], color = colors[2], width = 0.2, label='P3')
ax.bar(X + 0.6, iters[3], color = colors[3], width = 0.2, label='P4')

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.xlabel('N')
plt.xticks(np.arange(0.3,3.31,1),variables)
plt.ylabel('Iterations')
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)
# plt.yticks(np.arange(0,701,100))
plt.ylim([0,350])
plt.rcParams.update({'font.size' : 14})
plt.legend(fontsize=14)
plt.savefig('iter_sin2.png',dpi=600,bbox_inches='tight')
plt.show()



# read HTN
f = open('out2b')
N = Ntrials
variables = []
errors = np.zeros((4,N),dtype=float)
times = np.zeros((4,N),dtype=float)
iters = np.zeros((4,N),dtype=int)
for i in range(N):
   for j in range(Nheader):
      variables.append(f.readline().split()[-1]) # skip the header
   for j in range(4):
      f.readline() # Problem
      f.readline() # Params
      errors[j][i] = f.readline().split()[-1] # Error
      times[j][i] = f.readline().split()[-1]  # WTime
      iters[j][i] = f.readline().split()[-1]  # Iters
f.close()

## runtime for HTN
fig = plt.figure(3, figsize=(6,6))
plt.plot(variables,times[0],lines[0],label='P1',marker=marks[0],mfc='w',color=colors[0],ms=8)
plt.plot(variables,times[1],lines[1],label='P2',marker=marks[1],mfc='w',color=colors[1],ms=8)
plt.plot(variables,times[2],lines[2],label='P3',marker=marks[2],mfc='w',color=colors[2],ms=8)
plt.plot(variables,times[3],lines[3],label='P4',marker=marks[3],mfc='w',color=colors[3],ms=8)

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)

plt.xlabel('N')
plt.ylabel('Runtime (sec)')
# plt.yticks(np.arange(0,111,10))
plt.ylim([0,225])
plt.xticks(variables)
plt.rcParams.update({'font.size' : 14})
plt.legend(fontsize=14)
plt.savefig('time_htn2.png',dpi=600,bbox_inches='tight')
plt.show()


## Iteartions for SIN
fig = plt.figure(4, figsize=(6,6))
ax2 = fig.add_axes([0,0,1,1])
X = np.arange(N)
ax2.bar(X + 0.0, iters[0], color = colors[0], width = 0.2, label='P1')
ax2.bar(X + 0.2, iters[1], color = colors[1], width = 0.2, label='P2')
ax2.bar(X + 0.4, iters[2], color = colors[2], width = 0.2, label='P3')
ax2.bar(X + 0.6, iters[3], color = colors[3], width = 0.2, label='P4')

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(axis='x', which='minor', bottom=False, top=False)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)


plt.xlabel('N')
plt.xticks(np.arange(0.3,3.31,1),variables)
plt.ylabel('Iterations')
# plt.yticks(np.arange(0,701,100))
plt.ylim([0,350])
plt.rcParams.update({'font.size' : 14})
plt.legend(fontsize=14)
plt.savefig('iter_htn2.png',dpi=600,bbox_inches='tight')
plt.show()
