import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import subprocess
import argparse
import math
import sys

file = 'op.out'
n = 32

f = open(file)
errors = np.zeros((n),dtype=float)
times  = np.zeros((n),dtype=float)
iters  = np.zeros((n),dtype=int)

for i in range(n):
   f.readline() # op
   f.readline() # converged reason
   f.readline() # Problem
   f.readline() # Params
   errors[i] = f.readline().split()[-1] # Error
   times[i] = f.readline().split()[-1]  # WTime
   iters[i] = f.readline().split()[-1]  # Iters
f.close()

# print table
for i in range(n):
   print (" {:2d}  {:6.3e}   {:4d}".format(i+1,errors[i],iters[i]) )