include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

test1: test1.o MAfunctions.o
	-${CLINKER} -o test1 test1.o MAfunctions.o ${PETSC_LIB}
	${RM} test1.o MAfunctions.o


bratu: bratu.c 
	/usr/bin/mpicc $^ -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3 -O0  -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3 -O0    -I/home/takahashi/petsc/include -I/home/takahashi/petsc/linux-gnu-c-debug/include -I/usr/include -Wl,-rpath,/home/takahashi/petsc/linux-gnu-c-debug/lib -L/home/takahashi/petsc/linux-gnu-c-debug/lib -Wl,-rpath,/opt/intel/oneapi/mkl/latest/lib/intel64 -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,-rpath,/opt/intel/oneapi/compiler/2021.1.2/linux/lib -L/opt/intel/oneapi/compiler/2021.1.2/linux/lib -Wl,-rpath,/opt/intel/oneapi/mpi/2021.1.1/libfabric/lib -L/opt/intel/oneapi/mpi/2021.1.1/libfabric/lib -Wl,-rpath,/opt/intel/oneapi/mpi/2021.1.1/lib -L/opt/intel/oneapi/mpi/2021.1.1/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64_lin -Wl,-rpath,/opt/intel/oneapi/mpi/2021.1.1/lib/release -L/opt/intel/oneapi/mpi/2021.1.1/lib/release -Wl,-rpath,/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -L/opt/intel/oneapi/tbb/2021.1.1/lib/intel64/gcc4.8 -lpetsc -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -lstdc++ -ldl -lmpichfort -lmpich -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl -o $@


# etc
.PHONY: distclean test1

distclean:
	@rm -f *~ test1 *tmp