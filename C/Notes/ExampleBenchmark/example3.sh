#!/bin/bash

# run me with `bash example3.sh`
# make sure you first do `chmod +x example3.sh``

rm -f ex3.out

for i in {1..8}
do
   echo "n = $i" >> "ex3.out"
   /usr/bin/time -o "ex3.out" -a -p mpiexec -n $i ./tri -tri_m 1000000
done

# When I run this, my best time is n = 2: real 1.32, user 2.27, sys 0.30.