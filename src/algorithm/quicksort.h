#include "../datastruct/includeds.h"
#include <stdio.h>
//Implementation for points
//Swaps elements
static void swapPoint(double *a,double *b){
    //printf("%lf\n",a->x);
    double temp = *a;
    *a = *b;
    *b = temp;
}

static int partition(struct point **points,int dim,int left,int right){
    printf("l and r :%d %d\n",left,right);
    //Determine Axis To sort
    double *pivot,*comp;
    if(dim == 0){
        double *pivot = points[right]->x;
        double *comp = points[j]->x;
    }
    else{
        if(dim == 1){
            double pivot = points[right]->x;
            double comp = points[j]->x;
        }
        else{
            double pivot = points[right]->x;
            double comp = points[j]->x;
        }
    }
    int j;
    int i = left-1;
    for(j = left; j < right;j++){    
    //Gets compare axis
        if(comp < pivot){
            printf("compx,pivx: %lf %lf\n",comp,pivot);
            i++;
            printf("i and j :%d %d\n",i,j);
            struct point *swapi = points[i];
            swapPoint(*(swapi->x),points[j]->x);
            swapPoint(points[i],points[j]);
            swapPoint(points[i],points[j]);
        }
    }
    swapPoint(points[i+1],&points[right]);
    return (i+1);
}
void quicksortPoint(struct point **points,int dim,int left,int right){
    if (left < right){
        int pi = partition(points,dim,left,right);
        quicksortPoint(points,dim,left,pi-1);
        quicksortPoint(points,dim,pi+1,right);
    }
}
