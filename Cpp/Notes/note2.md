# Note 2

This is about running jobs on the cluster. Check out the subfolder `Lochness`.


## Getting started
You can access the cluster through Visual Studio Code with an SSH extension. Intall VS Code and then install the official SSH extension by Microsoft. Then hit `F1` and type `ssh` and scroll down to the first option: `Remote-ssh: connect to host...`. The address that you want to connect to is `ssh ucid@lochness.njit.edu`.

To connect the Lochness cluster, you have to on the NJIT network. You can use the NJIT Cisco VPN if you are off campus.

1. Type `ssh ucid@lochness.njit.edu` and then enter your UCID password when prompted.
2. Then type in `git pull https://github.com/tt73/MA-DDM.git` to wherever you want this repo downloaded.

## Modules
The cluster has many libraries installed but you have to first load them. To see what modules are available type `module avail`. You will not see PETSc on that list. To search for specific ones type `module spider petsc`. It will reveal that it actually has `petsc/3.12.0` available. Type in `module spider petsc/3.12.0` to see how the module can be loaded.

Usually you just type `module load python` to load up a module but petsc has many prerequisites:

   gnu8/8.3.0  impi/2019.4.243
   gnu8/8.3.0  impi/2019.6.166
   gnu8/8.3.0  mpich/3.3.1
   gnu8/8.3.0  mvapich2/2.3.2
   gnu8/8.3.0  openmpi3/3.1.4
   intel/19.0.4.243  impi/2019.4.243
   intel/19.0.4.243  impi/2019.6.166
   intel/19.0.4.243  mpich/3.3.1
   intel/19.0.4.243  mvapich2/2.3.2
   intel/19.1.0.166  impi/2019.4.243
   intel/19.1.0.166  impi/2019.6.166
   intel/19.1.0.166  mpich/3.3.1
   intel/19.1.0.166  mvapich2/2.3.2

This output makes no sense. You only need one of these combos:
* gnu8, mpich
* intel mpich
You can't load gnu8 and intel at the same time.

## Running jobs
You need to write a shell script and then submit the script to a cluster manager. The script `Lochness/lochjob` is configured to run the python script on the cluster. You can submit the script by typing `sbatch lochjob`. When the job is done, you will get an email.


