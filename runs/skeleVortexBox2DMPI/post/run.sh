#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

if [ $# -eq 0 ]
  then
    time_start=0 #t = 0.000
    time_end=5000 #t = 5.0
else
    time_start="$1"
    time_end="$2"
fi

clear
./clean.sh

./compile.sh
inc=250
D=$time_start
while [ $D -le $time_end ]; do
    newend=$((D + inc))
    echo $D $newend
    mpirun -np 2 ./drop 9 $D $newend
    D=$((D + inc))
done
echo "finished"
#run in mpi
#mpirun -np 2 ./drop
#mpirun -np 2 xterm -e gdb -ex run --args ./drop 9 time_start time_end
