#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

CC99='mpicc -std=c99' qcc -Wall -grid=quadtree -O2 -D_MPI=1 drop.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL  -lm





