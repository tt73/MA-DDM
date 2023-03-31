
## Let's try running the HTN DDM with a global tol of 1e-8 and 1e-4.
## We compare the error from both methods.


N=200
nd=4
ntol=1e-3
ktol=1e-2
op=0.2

gtol=1e-4

mpiexec -n $nd  ../maddm -N $N -htn -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -snes_rtol $gtol -sol
mv load_exact.m load_exact4.m
mv load_u.m load_u4.m

gtol=1e-8

mpiexec -n $nd  ../maddm -N $N -htn -sub_snes_rtol $ntol -sub_ksp_rtol $ktol -op $op -snes_rtol $gtol -sol
mv load_exact.m load_exact8.m
mv load_u.m load_u8.m