#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

#Normal Run
#np=0
#python3 plot.py $np 2 0 100 0
#./movie.sh


#MPI Run
np=0
python3 plot.py $np 2 0 100 1
./movie.sh



