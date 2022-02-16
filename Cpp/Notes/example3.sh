#!/bin/bash

# run me with `sh example3.sh`
# make sure you first do `chmod +x example3.sh``
echo "result for n = 1, m = 1000"
time mpiexec -n 1 ./tri -tri_m 1000

for i in {1..8..2}
do
   echo "n = $i"
   # time mpiexec -n $i ./tri -tri_m 100000
done