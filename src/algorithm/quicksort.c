#include "quicksort.h"
#include <stdio.h>
//Implementation for points
//Swaps elements
void swapPoint(struct point *a,struct point *b){
    struct point temp = *a;
    *a = *b;
    *b = temp;
}

int partition(struct point *points[],int *dim,int *left,int *right){
    //Determine Axis To sort
    double pivot;
    if(*dim == 0){
        pivot = points[*right]->x;
    }
    else{
        if(*dim == 1)
            pivot = points[*right]->y;
        else{
            pivot = points[*right]->z;
        }
    }
    int j;
    int i = *left - 1; 
    for(j = *left; j < *right;j++){    
        //Gets compare axis
        double comp;
        if(*dim == 0){
            comp = points[j]->x;
        }
        else {
            if(*dim == 1){
                comp = points[j]->y;
            }
            else{
                comp = points[j]->z;
            }
        }
        if(comp <= pivot){
            i++;
            swapPoint(points[i],points[j]);
        }
    }
    swapPoint(points[i+1],points[*right]);
    return (i+1);
}
void quicksortPoint(struct point *points[],int *dim,int *left,int *right){
    printf("going\n");
    if (left < right){
        printf("sorting...\n");
        //Partition list
        int pi = partition(points,dim,left,right);
        quicksortPoint(points,dim,left,&pi-1);
        quicksortPoint(points,dim,&pi+1,right);
    }
    else{
        printf("error\n");
    }
}
