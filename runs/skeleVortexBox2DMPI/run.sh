#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#serial


#./vortex 10 650

#run in mpi
mpirun -np 12 ./vortex 11 650
#gdb
#mpirun -np 12 xterm -e gdb -ex run --args ./vortex 11 650
#valgrind
#mpirun -np 12 valgrind --tool=memcheck --leak-check=full --log-file=valgrind.txt ./vortex 10 650
