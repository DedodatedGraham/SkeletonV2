#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#run in mpi
#mpirun -np 2 ./drop
mpirun -np 2 xterm -e gdb -ex run --args ./drop 9 650
