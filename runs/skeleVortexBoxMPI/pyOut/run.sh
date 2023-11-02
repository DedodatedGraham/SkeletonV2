#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

#Normal Run
#np=4
#python3 plot.py $np 2 0 15000 0
#./movie.sh


#MPI Run
np=4
python3 plot.py $np 2 0 15000 1
./movie.sh



