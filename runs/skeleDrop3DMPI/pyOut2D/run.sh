#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

#Normal Run
#np=8
#python3 plot.py $np 3 0 2 0
#./movie.sh


#MPI Run
np=3000
python3 plot.py $np 3 0 1 1
./movie.sh



