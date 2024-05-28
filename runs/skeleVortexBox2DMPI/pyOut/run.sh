#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

./clean.sh

#Normal Run
#np=1
#python3 plot.py $np 2 0 1000 0
#./movie.sh


#MPI Run
np=12
python3 plot.py $np 2 0 1800 1
./movie.sh $1
