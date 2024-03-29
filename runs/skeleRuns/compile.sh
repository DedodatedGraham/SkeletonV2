#rm _drop.c drop
#serial
#qcc -Wall -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -O2 drop.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
#test
#qcc -g -Wall -grid=quadtree -std=c99 -DFILTERED=1 -DMTRACE=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

#qcc -Wall -O2 drop.c -o drop -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm


#qcc -g -Wall -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL -lm

#openmp
qcc -g -Wall -fopenmp -grid=quadtree -std=c99 -DFILTERED=1 -D_UNREFINE=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -I$GSL -lgsl -lgslcblas -L$GSL -lm

#mpi
#qcc -source -grid=quadtree -D_MPI=1 drop.c
#mpicc -Wall -std=c99 -O2 _drop.c -o drop_mpi -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
