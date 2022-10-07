import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import subprocess
import argparse
import math
import sys


text = 'This script plots 1D solutions and saves it into the workspace.'
parser = argparse.ArgumentParser(description=text)
parser.add_argument("-n", "--num_nodes", help="number of nodes in the solution",type=int)
parser.add_argument("-p", "--problem", help="problem number 1, 2, 3, or 4",type=int)
args = parser.parse_args()

if args.num_nodes:
   n = args.num_nodes
else:
   n = 4
if args.problem:
   p = args.problem
else:
   p = 1

if p == 1:
   h = 2.0/(n+1)
   x = np.linspace(-1+h,1-h,n)
elif p == 2:
   h = 1.0/(n+1)
   x = np.linspace(h,1-h,n)
elif p == 3:
   h = 1.0/(n+1)
   x = np.linspace(h,1-h,n)
elif p ==4:
   h = 2.0/(n+1)
   x = np.linspace(-1+h,1-h,n)
else:
   print("invalid p")
   exit

# read the solution
u = np.zeros(n)
f = open("load_u.m")
f.readline()
f.readline()
f.readline()
for i in range(n):
   u[i] = f.readline()
f.close()

# read the exact solution
u_exact = np.zeros(n)
f = open("load_exact.m")
f.readline()
f.readline()
f.readline()
for i in range(n):
   u_exact[i] = f.readline()
f.close()


# plot settings
plt.rcParams.update({'font.size' : 14})
colors=sns.color_palette("rocket",2)
fig = plt.figure(1, figsize=(6,6))

plt.plot(x,u,'-', label='numerical',marker='.',mfc='w',color=colors[0],ms=8)
# plt.plot(x,u,'-', label='numerical',color=colors[0],ms=8)
plt.plot(x,u_exact,'--', label='exact',color=colors[1],ms=8)


plt.tick_params(direction='in',right=True,top=True)
plt.tick_params(labelsize=14)
plt.tick_params(direction='in',which='minor', length=5, bottom=True, top=True, left=True, right=True)
plt.tick_params(direction='in',which='major', length=10, bottom=True, top=True, left=True, right=True)
plt.tick_params(labelbottom=True,labeltop=False,labelright=False,labelleft=True)
plt.legend(fontsize=14)
plt.xlabel('x')
plt.ylabel('u(x)')
plt.savefig('sol.png',dpi=600,bbox_inches='tight')
plt.show()

