#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

#Normal Run
#np=12
#python3 plot.py $np 3 0 250 0
#./movie.sh


#MPI Run
np=12
python3 plot.py $np 3 0 250 1
./movie.sh



