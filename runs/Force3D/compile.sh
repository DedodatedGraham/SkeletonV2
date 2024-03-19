#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


#compile serial
qcc -g -Wall -grid=octree -O2 -DMTRACE=0 cyl.c -o cyl -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

#compile openmp
#qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -DMTRACE=0 -O2 drop-Skele.c -o _drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL -lm

#compile MPI
#qcc -source -grid=octree -D_MPI=1 -DMTRACE=0 cyl.c
#mpicc -g -Wall -std=c99 -D_XOPEN_SOURCE=700 -O2 _cyl.c -o cyl -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
