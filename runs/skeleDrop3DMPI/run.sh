#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#run in mpi

#gdb
mpirun -np 2 xterm -e gdb -ex run --args ./vortex 7 650
#valgrind
#mpirun -np 2 xterm -e valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=valgrind.txt ./vortex 7 650
