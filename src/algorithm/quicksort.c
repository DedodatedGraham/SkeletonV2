#include "quicksort.h"
#include <stdio.h>
//Implementation for points
//Swaps elements
void swapPoint(struct point *a,struct point *b){
    //printf("%lf\n",a->x);
    struct point *temp = a;
    a = b;
    b = temp;
}

int partition(struct point *points[],int *dim,int left,int right){
    //Determine Axis To sort
    struct point *pivot = points[right];
    int j;
    int i = left - 1;
    for(j = 0; j < right;j++){    
        //Gets compare axis
        struct point *comp = points[j];
        if(*dim == 0){
            if(&comp->x <= &pivot->x){
                i++;
                swapPoint(points[i],points[j]);
            }
        }
        else {
            if(*dim == 1){
                if(&comp->y <= &pivot->y){
                    i++;
                    swapPoint(points[i],points[j]);
                }
            }
            else{
                if(&comp->z <= &pivot->z){
                    i++;
                    swapPoint(points[i],points[j]);
                }
            }
        }
    }
    swapPoint(points[i+1],points[right]);
    return (i+1);
}
struct point *quicksortPoint(struct point points[],int *dim,int left,int right){
    if (left < right){
        //Partition list
        int pi = partition(&points,dim,left,right);
        struct point *lefts;
        struct point *rights;
        lefts = quicksortPoint(points,dim,left,pi-1);
        rights = quicksortPoint(points,dim,pi+1,right);
        int i;
        struct point *returns[right-left];
        for(i = 0; i < (right - left);i++){
            if(i < pi-1){
                *returns[i] = lefts[i];
            }
            else if (i > pi){
                *returns[i] =  rights[i - (pi)];
            }
            else{
                *returns[i] = points[pi];
            }
        }
        return *returns;
    }
    return NULL;
}
