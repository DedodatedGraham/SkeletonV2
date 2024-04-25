#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"


#compile serial
#qcc -dimensions -g -Wall -grid=quadtree -O2 -DMTRACE=0 vortex.c -o vortex -I$BASILISK/gl -L$BASILISK/gl -lfb_tiny -lglutils -lm
#qcc -g -Wall -grid=quadtree -O2 vortex.c -o vortex -I$BASILISK/gl -L/$BASILISK/gl -lm

#-lOSMesa -lGLU -lOSMesa
#compile openmp
#qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -DMTRACE=0 -O2 drop-Skele.c -o _drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL -lm

#compile MPI
CC99='mpicc -std=c99' qcc -g -Wall -grid=quadtree -O2 -D_MPI=1 -DMTRACE=1 vortex.c -o vortex -lm





