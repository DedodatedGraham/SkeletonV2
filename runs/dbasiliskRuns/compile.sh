#rm _drop.c drop

qcc -g -Wall -grid=quadtree -std=c99 -DFILTERED=1 -DMTRACE=0 -D_UNREFINE=1 -O1 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$BASILISK/skele/gsl-2.7.1 -lgsl -lgslcblas -L$BASILISK/skele/gsl-2.7.1 -lm
#qcc -g -Wall -grid=quadtree -std=c99 -DFILTERED=1 -DMTRACE=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

#openmp
#qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

#mpi
