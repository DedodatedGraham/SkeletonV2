#include "quicksort.h"
#include <stdio.h>
//Implementation for points
//Swaps elements
void swapPoint(struct point *a,struct point *b){
    //printf("%lf\n",a->x);
    struct point temp = *a;
    *a = *b;
    *b = temp;
}

int partition(struct point **points[],int *dim,int left,int right){
    printf("l and r :%d %d\n",left,right);
    //Determine Axis To sort
    struct point pivot = **points[right];
    int j;
    int i = left-1;
    for(j = left; j < right;j++){    
    //Gets compare axis
    struct point comp = **points[j];
        if(*dim == 0){
            if(comp.x < pivot.x){
                printf("compx,pivx: %lf %lf\n",comp.x,pivot.x);
                i++;
                printf("i and j :%d %d\n",i,j);
                swapPoint(*points[i],*points[j]);
            }
        }
        else {
            if(*dim == 1){
                if(comp.y <= pivot.y){
                    i++;
                    printf("i and j :%d %d\n",i,j);
                    swapPoint(*points[i],*points[j]);
                }
            }
            else{
                if(comp.z <= pivot.z){
                    i++;
                    printf("i and j :%d %d\n",i,j);
                    swapPoint(*points[i],*points[j]);
                }
            }
        }
    }
    swapPoint(*points[i+1],*points[right]);
    return (i+1);
}
void quicksortPoint(struct point **points[],int *dim,int left,int right){
    if (left < right){
        //Partition list
        printf("Finding Partition\n");
        int pi = partition(points,dim,left,right);
        printf("Found %d\n",pi);
        quicksortPoint(points,dim,left,pi-1);
        quicksortPoint(points,dim,pi+1,right);
    }
}
