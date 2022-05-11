# Note 1

## Goal
The goal is to be able to call to construct and solves systems all in petsc on the backend, and input problem parameters and visualize plots with Matlab on the front end.

## Install PETSc
Getting PETSc is actually pretty easy, especially on Linux. Go to this [page](https://petsc.org/release/download/) and follow instructions. Here are the barebones instructions:


1. Install MPICH on your computer. On Ubuntu, the command was `sudo apt install mpich`. It should be saved in `/usr/bin/mpich` by default. 
   
2. Optional: you can install the Intel Math Kernel Library. The BLAS/LAPACK library is a prereq for PETSc but it will MKL.  
   
3. Open up the terminal `Ctrl+Shift+t`. Type `cd` to move to your home directory. (or wherever you want). Then type `git clone -b release https://gitlab.com/petsc/petsc.git petsc` to pull their source code.

4. Type `cd petsc` to move into the directory. And then type `./configure` to set up the code. This will work if you have MPI and BLAS/LAPACK istalled. You have the option to let PETSc install the dependencies for you. Do `./configure --with-mpi-dir=/usr --download-fblaslapack` and it should work. If you have MKL installed then omit the last argument, ie. `./configure --with-mpi-dir=/usr`. 

5. When the configure is done, scroll up in the termial a bit and you should see a message telling you to run a `make` command. Copy & paste that command into the terminal to install PETSc.

The installation is done.

## Additional installation steps

You need to be able to link the PETSc library when compiling code. To avoid manually writing gigantic compile commands, PETSc devs suggest you use Makefiles. They have written variables and compile command templates for you to import into your Makefile. But to find these files, you still need to tell the computer where PETSc is installed. To do this, you must create environment variables.

* In the terminal type `cd; gedit .bashrc` to open up an editor. Anywhere in the file, add the lines below and save the file. You have to modify the first line.
   ```
   export PETSC_DIR=/wherever/you/installed/petsc
   export PETSC_ARCH=linux-gnu-c-debug
   ```


## Compiling the code
PETSc devs want you to compile code with Makefiles. Create a file named `makefile` in the directory containing your source code. In your Makefile, add the lines:
   ```
   include ${PETSC_DIR}/lib/petsc/conf/variables
   include ${PETSC_DIR}/lib/petsc/conf/rules
   ```
These lines will import variables and compilation rules provided by the PETSc devs. They have already written rules to compile Fortran, C, and C++ code. All you have to do is type into the command line `make` followed by the name of your source code. For example, `make ex50.c` will compile the code and produce `ex50`.

If you have a code that has multiple dependencies e.g. a code that depends on several custom libraries, then you acutally need to write a custom rule to get the code compiling.

## Running the code
You run PETSc code with `mpirun`, `mpiexec`, or `mpifort` depending on your MPI implementation and source code. For this project, we are using C code so `mpiexec` is what we use if we want to call it from the command line. 
The order of the arguments matter. Typical command or running PETSc code might look something like this:
```
   mpiexec -n 4 ./program_name -flag_for_program input
```
Notice that the `mpiexec` command and the program you wish to run both take inputs. Inputs for the program come after the name of the program. A common mpiexec argument is the number of threads `-n`. This has to come before the program name. The program has to have the `./` to signify that the program you are trying to call lives in the current directory. The program might also have a flag `-n` which can be confusing.