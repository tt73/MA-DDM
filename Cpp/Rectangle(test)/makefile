include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

test1: test1.o MAfunctions.o
	-${CLINKER} -o test1 test1.o MAfunctions.o ${PETSC_LIB}
	${RM} test1.o MAfunctions.o

# etc
.PHONY: distclean test1

distclean:
	@rm -f *~ test1 *tmp