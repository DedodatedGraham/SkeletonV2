#serial
#qcc -g -Wall -grid=quadtree -std=c99 -DFILTERED=1 -DMTRACE=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_dumb -lm

#OpenMP
qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -DMTRACE=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_dumb -lm

#OpenMPI
#qcc -source -grid=quadtree -D_MPI=1 dropslide.c
#mpicc -Wall -std=c99 -O2 _dropslide.c -o drop_mpi -I$BASILISK -L$BASILISK/gl -lglutils -lfb_dumb -lm





