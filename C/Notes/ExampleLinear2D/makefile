# -*- mode: makefile -*-

#  This sample Makefile can be used to compile PETSc applications
#  Copy this file to your source directory as "Makefile" and modify as needed.
#  See also $PETSC_DIR/share/petsc/Makefile.user for the preferred approach
#  You must set the environmental variable(s) PETSC_DIR (and PETSC_ARCH if PETSc was not configured with the --prefix option)
#
#   You can set specific values below but that should rarely need to
CFLAGS		 =
FFLAGS		 =
CPPFLAGS         =
FPPFLAGS         =

#  For example - a single source file (ex1.c or ex1.F90) can be compiled with:
#
#      make ex1
#
#  You do not need to edit this makefile at all.
#
#  For a multi-file case, suppose you have the source files a.c, b.c, and c.cxx
#  This can be built by uncommenting the following two lines.
#
# app : a.o b.o c.o
# 	$(LINK.C) -o $@ $^ $(LDLIBS)
#
#  When linking in a multi-files with Fortran source files a.F90, b.c, and c.cxx
#  You may need to use
#
# app : a.o b.o c.o
# 	$(LINK.F) -o $@ $^ $(LDLIBS)

# If the file c.cxx needs to link with a C++ standard library -lstdc++ , then
# you'll need to add it explicitly.  It can go in the rule above or be added to
# a target-specific variable by uncommenting the line below.
# app : LDLIBS += -lstdc++

include ${PETSC_DIR}/lib/petsc/conf/variables

#  To access the PETSc variables for the build, including compilers, compiler flags, libraries etc but
#  manage the build rules yourself (rarely needed) comment out the next lines
#
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test


# Type `make all` the build the two programs in this directory
all: poisson view2d

# Type `make clr` to delete the executable files.
clr:
	rm -f poisson view2d
	rm -f load_b.m load_u.m load_mat.m view2d.out