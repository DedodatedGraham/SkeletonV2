#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#serial


#./vortex 10 650

#run in mpi
#gdb
mpirun -np 2 xterm -e gdb -ex run --args ./vortex 10 650
#valgrind
#mpirun -np 2 xterm -e valgrind --tool=memcheck --leak-check=full --log-file=valgrind.txt ./vortex 10 650
