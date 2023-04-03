#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
./clean.sh
./compile.sh
export OMP_NUM_THREADS=12
./drop
#gdb -q -ex=r drop
#./../pScript/run.sh
