## uncomment below if compiling on cluster
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
# CFLAGS += -pedantic -std=c99

maddm: maddm.o MAfunctions.o
	-${CLINKER} -o maddm maddm.o MAfunctions.o ${PETSC_LIB}
	${RM} maddm.o MAfunctions.o

bratu: bratu.c
	${CC} $^ ${PETSC_LIB} ${PETSC_CCPPFLAGS} ${CC_FLAGS} -o $@

pbratu: pbratu.c
	${CC} $^ ${PETSC_LIB} ${PETSC_CCPPFLAGS} ${CC_FLAGS} -o $@


# etc
.PHONY: distclean test1

distclean:
	@rm -f *~ test1 *tmp
