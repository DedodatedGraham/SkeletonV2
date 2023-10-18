#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

clear
./clean.sh

np=2
python3 plot.py $np 2 651 800
./movie.sh



