# Note 1 

## Goal 
The goal is to be able to call to construct and solves systems all in petsc on the backend, and input problem parameters and visualize plots with Matlab on the front end. 

## Install PETSc 
Getting PETSc is actually pretty easy, especially on Linux. Go to this [page](https://petsc.org/release/download/) and follow instructions. Here are the barebones instructions: 
1. Open up the terminal `Ctrl + Shift + t` 
2. Type `cd` to move to your home directory. (or wherever you want)
3. Type `git clone -b release https://gitlab.com/petsc/petsc.git petsc` to pull their source code 
4. Type `cd petsc` to move into the directory. And then type `./configure` to set up the code. It may take a while. If it fails, then you may be lacking prerequisites. If its a basic prereq like MIPCH then the installer already has a workaround. Do `./configure --download-mpich` and it should work. 
5. After the `./configure` is done, scroll up in the termial a bit and you should see some instructions in the output. Run the provided `make` command to install PETSc. 
6. The configure file will tell what your environments variables are needed. You need to have variables `PETSC_DIR` and `PETSC_ARCH` defined in the shell. In the terminal type `cd; gedit .bashrc` to open up an editor. Add the lines below and save the file. 
   ```
   export PETSC_DIR=/wherever/you/installed/petsc
   export PETSC_ARCH=linux-gnu-c-debug
   ```

## Compiling the code 

