## First load the python module with `module load python`.
## Then gnenerate the job scripts with `python generate_tests.py`
## Then subit the jobs on the cluster with this scrpit with `bash submit_tests`

# delete all out files
rm out*

# submit all the jobs
for N in 2 4 6 8 9
do
   printf "Running N = $N... \n"
   bash job$N
done
