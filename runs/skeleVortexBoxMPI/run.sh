#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#run in mpi
mpirun -np 4 xterm -e gdb -ex run --args ./vortex 7 650
