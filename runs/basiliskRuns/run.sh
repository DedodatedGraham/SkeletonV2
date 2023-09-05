#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"
clear
if [ $# -eq 0 ]
  then
    time_start=0
    threads=14
else
  if [ $# -eq 1 ]
    then
      time_start=$1
      threads=14
  else
    time_start=$1
    threads=$2
  fi
fi
if [ $time_start -eq 0 ]
  then
    ./clean.sh
    echo "cleaned"
fi
./compile.sh
export OMP_NUM_THREADS=$threads
./drop $time_start
#gdb -q -ex=r drop
#./../pScript/run.sh
