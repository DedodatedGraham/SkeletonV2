#!/bin/bash

parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
cd "$parent_path"

qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -O2 drop-NoSkele.c -o _drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL -lm
