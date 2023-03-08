import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# data settings
Ntrials = 8
Nvariables = 1
Nheader = 1
fname  = 'out_N'
var_labels = ["HTN","SIN","NKS"]

f = open(fname)
N = Ntrials
variables = []
errors = np.zeros((Nvariables,Ntrials),dtype=float)
times = np.zeros((Nvariables,Ntrials),dtype=float)
iters = np.zeros((Nvariables,Ntrials),dtype=int)

# Loop over all trials
#    Skip the header lines
#    Loop over the 4 examples
for i in range(Ntrials):
   for j in range(Nheader):
      variables.append(str(f.readline().split()[-1])) # skip the header
   for j in range(Nvariables):
      f.readline() # Problem
      f.readline() # Params
      errors[j][i] = f.readline().split()[-1] # Error
      times[j][i] = f.readline().split()[-1]  # WTime
      iters[j][i] = f.readline().split()[-1]  # Iters
f.close()

## print out tables

print("runtime")
rowstr = ""
for i in range(Ntrials):
   rowstr += "{:>12s} &".format(variables[i])
print(rowstr)
for j in range(Nvariables):
   rowstr = ""
   for i in range(Ntrials):
      rowstr += "{:12.3f} &".format(times[j,i])
   print(rowstr)

print("error")
rowstr = ""
for i in range(Ntrials):
   rowstr += "{:>12s} &".format(variables[i])
print(rowstr)

for j in range(Nvariables):
   rowstr = ""
   for i in range(Ntrials):
      rowstr += "{:12.3E} &".format(errors[j,i])
   print(rowstr)

print("iterations")
rowstr = ""
for i in range(Ntrials):
   rowstr += "{:>12s} &".format(variables[i])
print(rowstr)
for j in range(Nvariables):
   rowstr = ""
   for i in range(Ntrials):
      rowstr += "{:12d} &".format(iters[j,i])
   print(rowstr)



plt.rcParams.update({'font.size' : 14})
colors = sns.color_palette("hls",Nvariables)
lines = ['-',':','--','.-']
marks = ['o','s','^','x']

### DO LINE PLOT FOR THE RUNTIME
fig = plt.figure(1, figsize=(6,6))
for i in range(Nvariables):
   plt.plot(variables,times[i],lines[i],label=var_labels[i],marker=marks[i],mfc='w',color=colors[i],ms=8)

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)

plt.xlabel('Overlap Percentage')
plt.ylabel('Runtime (sec)')
# plt.yticks(np.arange(0,111,10))
plt.xticks(variables)

plt.legend(fontsize=14)
plt.savefig('runtime.png',dpi=200,bbox_inches='tight')
plt.show()

## DO BAR GRAPH FOR THE ITERATIONS
fig = plt.figure(2, figsize=(6,6))
ax = fig.add_axes([0,0,1,1])
X = np.arange(N)
dx = 1.0/Nvariables
for i in range(Nvariables):
   ax.bar(X + i*dx, iters[i], color = colors[i], width = dx, label=var_labels[i])
# ax.bar(X + 0.2, iters[1], color = colors[1], width = 0.2, label='P2')
# ax.bar(X + 0.4, iters[2], color = colors[2], width = 0.2, label='P3')
# ax.bar(X + 0.6, iters[3], color = colors[3], width = 0.2, label='P4')

plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.xlabel('Overlap Percentage')
plt.xticks(np.arange(0.3,Ntrials+0.3,1),variables)
plt.ylabel('Iterations')
plt.minorticks_on(); plt.tick_params(axis='x', which='minor', bottom=False, top=False)
# plt.yticks(np.arange(0,701,100))
# plt.ylim([0,650])
plt.rcParams.update({'font.size' : 14})
plt.legend(fontsize=14)
plt.savefig('iters.png',dpi=200,bbox_inches='tight')
plt.show()