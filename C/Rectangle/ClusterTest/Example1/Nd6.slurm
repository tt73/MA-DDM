Starting /home/t/tt73/.bash_profile ... standard AFS bash profile

========================
Home directory : /home/t/tt73 is not in AFS -- skipping quota check
========================

On host node004 :
	 09:12:34 up 333 days, 14:45,  0 users,  load average: 0.08, 0.03, 0.05

=== === === Your Kerberos ticket and AFS token status === === ===
klist: No credentials cache found (filename: /tmp/krb5cc_90502_XkEizAyZlO)
Kerberos : 
AFS      : 

Loading default modules ...
Create file : "/home/t/tt73/.modules" to customize.


Currently Loaded Modules:
  1) python/3.8   3) mpich/3.3.1    5) openblas/0.3.7    7) petsc/3.12.0
  2) gnu8/8.3.0   4) phdf5/1.10.5   6) scalapack/2.0.2

 

[-1.5, 1.5]^2
tol = 1e-6
Now running mpiexec -np 6 ../../maddm -N 60 -problem ex1 -sub_snes_rtol 1e-6 -sub_ksp_rtol 1e-6 -xmin -1.5 -xmax 1.5 -ymin -1.5 -ymax 1.5 -snes_converged_reason -snes_rtol 1e-9 
Now running mpiexec -np 6 ../../maddm -N 300 -problem ex1 -sub_snes_rtol 1e-6 -sub_ksp_rtol 1e-6 -xmin -1.5 -xmax 1.5 -ymin -1.5 -ymax 1.5 -snes_converged_reason -snes_rtol 1e-9 
[-2.0, 2.0]^2
tol = 1e-6
Now running mpiexec -np 6 ../../maddm -N 80 -problem ex1 -sub_snes_rtol 1e-6 -sub_ksp_rtol 1e-6 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-9 
Now running mpiexec -np 6 ../../maddm -N 400 -problem ex1 -sub_snes_rtol 1e-6 -sub_ksp_rtol 1e-6 -xmin -2.0 -xmax 2.0 -ymin -2.0 -ymax 2.0 -snes_converged_reason -snes_rtol 1e-9 
done