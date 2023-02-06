import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# filename = 'out_test3'

# Reading files usually depends on 3 things:
# 1. Filename.
# 2. Number of trials. A trial consists of a run of each of the 3 problems.
# 3. Header size. This is the number of lines to skip for each trial. Usually 1.

print("Enter filename: ")
filename = input()
print("Enter number of trials: ")
Ntrials = int(input())
# print("Enter header size: ")
Nheader = 1

f = open(filename)

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

print(variables)
print("Errors:")
print(errors)
print("Times:")
print(times)
print("Iters:")
print(iters)


##  Make line graph for runtime
plt.rcParams.update({'font.size' : 14})
colors = sns.color_palette("rocket",4)
lines = ['-',':','--','.-']
marks = ['o','s','^','x']

fig = plt.figure(1, figsize=(6,6))
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

plt.xlabel('N')
# plt.xlabel('Overlap Percentage')
plt.ylabel('Runtime (sec)')
plt.xticks(np.arange(0.3,5.5,1),variables)
plt.legend(fontsize=14)
plt.savefig('runtime.png',dpi=300,bbox_inches='tight')
plt.show()


## Make bar graph for iters
fig = plt.figure(2, figsize=(6,6))
X = np.arange(N)
ax2 = fig.add_axes([0,0,1,1])
ax2.bar(X + 0.0, iters[0], color = colors[0], width = 0.2, label='P1')
ax2.bar(X + 0.2, iters[1], color = colors[1], width = 0.2, label='P2')
ax2.bar(X + 0.4, iters[2], color = colors[2], width = 0.2, label='P3')
ax2.bar(X + 0.6, iters[3], color = colors[3], width = 0.2, label='P4')

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)

plt.xlabel('N')
# plt.xlabel('Overlap Percentage')
plt.ylabel('Newton Iterations')
plt.xticks(np.arange(0.3,5.5,1),variables)
plt.legend(fontsize=14)
plt.savefig('iters.png',dpi=300,bbox_inches='tight')
plt.show()