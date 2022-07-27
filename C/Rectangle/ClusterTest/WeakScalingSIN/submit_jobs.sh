## First load the python module with `module load python`.
## Then gnenerate the job scripts with `python generate_tests.py`
## Then subit the jobs on the cluster with this scrpit with `bash submit_tests`

# delete all out files
rm -f out*

# submit all the jobs
for N in 1 4 9
do
   sbatch job$N
done

# check the queue
squeue