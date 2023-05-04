## To compile this on Stheno or Lochness, type `module load gnu8 mpich petsc/3.12.0`. Then type `make maddm`.
## You can run the compiled code serially with `./maddm`.
## You can use `mpirun` or `mpiexec` to run in parallel.
## The NJIT cluster offers a more recent version of PETSc, but it fails to run to some parallel examples.

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

maddm: maddm.o MAfunctions.o
	-${CLINKER} -o maddm maddm.o MAfunctions.o ${PETSC_LIB}
	${RM} maddm.o MAfunctions.o

interp: interp_example.o MAfunctions.o
	-${CLINKER} -o interp interp_example.o MAfunctions.o ${PETSC_LIB}
	${RM} interp_example.o MAfunctions.o

bratu: bratu.c
	${CC} $^ ${PETSC_LIB} ${PETSC_CCPPFLAGS} ${CC_FLAGS} -o $@

pbratu: pbratu.c
	${CC} $^ ${PETSC_LIB} ${PETSC_CCPPFLAGS} ${CC_FLAGS} -o $@


# etc
.PHONY: distclean test1

distclean:
	@rm -f *~ test1 *tmp
