#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

./compile.sh

#serial


./cyl
#./pyOut/run.sh
#run in mpi
#mpirun -np 2 ./cyl
