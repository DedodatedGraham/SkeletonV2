#qcc -Wall -fopenmp -grid=quadtree -std=c99 -D_OscilDrop=1 -O2 drop.c -o drop_op -I$HOME -L$HOME/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
#qcc -Wall -fopenmp -grid=quadtree -std=c99 -D_OscilDrop=1 -O2 drop.c -o drop_op -I$BASILISK -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
qcc -Wall -fopenmp -grid=quadtree -std=c99 -D_OscilDrop=1 -O2 dropslide.c -o drop -I$BASILISK -L$BASILISK/gl -lglutils -lfb_dumb -lm
