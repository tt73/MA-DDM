import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys


text = 'This is a script that runs a ddm over multiple values of k'
parser = argparse.ArgumentParser(description=text)
parser.add_argument("-v", "--verbose", help="show what commands are being run", action="store_true")
parser.add_argument("-p", "--plot", help="make plots of real and imag parts of solution", action="store_true")
parser.add_argument("-m", "--method", help="choose the method, -m number_from_1_to_8",type=int)
args = parser.parse_args()


f = open("out_test1")
N = 5
times = np.zeros((N,3),dtype=float)
for i in range(N):
   f.readline() # N
   for j in range(3):
      f.readline() # Problem
      f.readline() # Params
      f.readline() # Error
      times[i][j] = f.readline().split()[1] # WTime
      f.readline() # Iters

print(times)


ax1.plot(angles,reu,'-',label=str(waves[i]))
ax1.title.set_text('')
# ax2.plot(angles,imu,'-',label=str(waves[i]))
# ax2.title.set_text('Imaginary')
ax1.plot(angles,np.real(u),'k:',label='_nolegend_')
# ax2.plot(angles,np.imag(u),'k:',label='_nolegend_')
ax1.legend()
# ax2.legend()
plt.suptitle("Solution on exterior R = {}".format(r2))
plt.savefig('wavenumber.png')
plt.show()


# # Choose your method here...
# num = 7
# if args.method:
#    num = args.method # or at runtime
# methods = [
#    'abc_direct',   # 1
#    'abc_scalar',   # 2
#    'abc_oscalar',  # 3
#    'abc_lb',       # 4
#    'abc_pade',     # 5
#    'abc_opade',    # 6
#    'arc_scalar',   # 7
#    'arc_pade',     # 8
#    'ie_direct',     # 9
#    'arc_scalar -curve 2',   # 10
#    'arc_scalar -curve 3'    # 11
#           ]
# method = methods[num-1]
# tag = method[0:3]
# if(num==9):
#    tag = 'arc'


# # Choose your mesh
# file =  'circle8h0.100'
# # file =  'circle8h0.050'
# # file =  'circle8h0.075'

# # file =  'concen2h0.100'
# # file =  'concen2h0.075'

# # Choose your wave range
# waves = [1, 1.5, 2, 2.5, 3, 3.5, 4]

# N = len(waves)
# status    = np.zeros((N),dtype=np.int8)  # 1, 2, or 3
# num_iters = np.zeros((N),dtype=np.uint32)
# res_final = np.zeros((N))
# error_r   = np.zeros((N),dtype=float)
# error_i   = np.zeros((N),dtype=float)

# # if plotting...
# if args.plot:
#    ax1 = plt.subplot(121) # real
#    ax2 = plt.subplot(122) # imag
#    r1 = 1
#    r2 = 3
#    angle_inc = 0

# for i in range(N):
#    # cmd = "./{} -debug 0 -show 0 -file {} -k {:.2f} -scalar {:.6f},{:.6f}".format(method, file, waves[i],waves[i]*0.5,-waves[i])
#    # cmd = "./{} -debug 0 -show 0 -file {} -k {:.2f} -scalar {:.6f},{:.6f} -arcg".format(method, file, waves[i],0,-waves[i])
#    cmd = "./{} -debug 0 -show 0 -file {} -k {:.2f}".format(method, file, waves[i])
#    if(args.verbose):
#       print("Now running: "+cmd)
#    out = subprocess.check_output(cmd, shell=True, universal_newlines=True)
#    if(args.verbose):
#       print(out)
#    status[i] = out.split()[1]
#    num_iters[i] = out.split()[3]
#    res_final[i] = out.split()[5]
#    if(status[i]== 1 or status[i]==3):
#       if(tag == "arc"):
#          error_r[i], error_i[i] = em.compute_error_exterior(k=waves[i])
#       elif(tag == "abc"):
#          error_r[i], error_i[i] = em.compute_error_annulus(1,3,waves[i],180,'sol_ext.csv')
#    else:
#       print("aaa")
#       error_r[i] = np.nan
#       error_i[i] = np.nan
#    if args.plot:
#       angles, reu, imu = em.read_exterior_solution('sol_ext.csv')
#       ax1.plot(angles,reu,'-',label=str(waves[i]))
#       ax1.title.set_text('Real')
#       ax2.plot(angles,imu,'-',label=str(waves[i]))
#       ax2.title.set_text('Imaginary')
#       if(tag == "abc"):
#          u = em.compute_solution_annulus(r1,r2,waves[i],angle_inc*pi/180+pi,angles)
#       elif(tag=="arc"):
#          u = em.compute_solution_exterior(r1,r2,waves[i],angle_inc*pi/180,angles)
#       ax1.plot(angles,np.real(u),'k:',label='_nolegend_')
#       ax2.plot(angles,np.imag(u),'k:',label='_nolegend_')


# print("File:   {}".format(file))
# print("Method: {}".format(method))
# print("|  k  | iter| Code| Re(Error)| Im(Error)|")
# print("|----:|----:|----:| --------:| --------:|")
# for i in range(N):
#    print("|{:5.1f}|{:5d}|{:5d}|{:10.3f}|{:10.3f}|".format(waves[i],num_iters[i],status[i],error_r[i],error_i[i]))

# if args.plot:
#    ax1.legend()
#    ax2.legend()
#    plt.suptitle("Solution on exterior R = {}".format(r2))
#    plt.savefig('wavenumber.png')
#    plt.show()