#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

#qcc -g -Wall -fopenmp -grid=quadtree -O2 -DMTRACE=1 drop.c -o drop -L$BASILISKGL -lfb_tiny -lglutils -lm

#CC99='mpicc -std=c99' qcc -g -Wall -grid=octree -O2 -D_MPI=1 -DMTRACE=0 drop.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa
qcc -source -grid=quadtree -D_MPI=1 -DMTRACE=0 drop.c
mpicc -g -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _drop.c -o drop -L$BASILISKGL -lfb_tiny -lglutils -lm
