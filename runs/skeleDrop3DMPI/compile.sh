#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


#compile serial
qcc -g -Wall -grid=octree -O2 -DMTRACE=0 drop.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL  -lm

#compile MPI
#CC99='mpicc -std=c99' qcc -g -Wall -grid=octree -O2 -D_MPI=1 -DMTRACE=1 drop.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL  -lm





