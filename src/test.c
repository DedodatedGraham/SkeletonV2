#include "skeletize.h"

int main(int argc, char *argv[]){
    //Load file
    char *path;
    if(argc >= 2){
        path = argv[1];
    }
    FILE *file;
    file = fopen(path,"r");
    double **points;
    //load in points
    int index = 0;
    while(fscanf(file, "%d %d %d %d", &points[index][0],&points[index][1],&points[index][2],&points[index][3]) == 4){
        index++;
    }
    skeletize(points);
}
