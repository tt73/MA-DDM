#!/bin/sh 
# to run me, type `source load_petsc.sh`
module load  gnu8
module load  mpich 
module load petsc/3.12.0

# If petsc was sucessfully loaded then the following command will show petsc. 
module list

# See if PETSC_DIR and PETSC_ARCH are defined.  
echo "PETSC_DIR:"  $PETSC_DIR 
echo "PETSC_ARCH:" $PETSC_ARCH