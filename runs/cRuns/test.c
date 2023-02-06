#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../../src/skeletize.h"

int main(){
    
    //open file
    FILE *fp1 = fopen("../../indata/infc-0.000.dat","r");
    if(fp1 == NULL){
        fprintf(stdout,"error opening file\n");
        return 1;
    }
    int sz = 0;
    int ch = 0;
    double tx;
    double ty;
    double nx;
    double ny;
    //measure the length of input file
    while((ch = fgetc(fp1)) != EOF){
        if(ch == '\n'){
            sz ++;
        }
    } 

    //alllocate for load in properties
    rewind(fp1);
    double **inputpoints = (double**)malloc(sz*sizeof(double*));
    int count = 0;
    fprintf(stdout,"trying length %d\n",sz); 
    
    //scan file
    while(count < sz){
        fscanf(fp1,"%lf %lf %lf %lf",&tx,&ty,&nx,&ny);
        inputpoints[count] = (double*)malloc(4*sizeof(double));
        inputpoints[count][0] = tx;
        inputpoints[count][1] = ty;
        inputpoints[count][2] = nx;
        inputpoints[count][3] = ny;
        fprintf(stdout,"loaded in x:%f y:%f nx:%f ny:%f\n",tx,ty,nx,ny);
        count++;
    }

    //close
    fclose(fp1);
    int amount = 2;
    sz = sz - 1;
    char name[80];
    sprintf(name,"testOut.dat");
    skeletize(inputpoints,&sz,&amount,name);
}
