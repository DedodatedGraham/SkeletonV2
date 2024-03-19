#include <stdio.h>
#include <math.h>
#include <time.h>

void skeleSwap(double *a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int skelePartition(double **arr, int lbound, int tbound, int axis,int mode){
    double pivot = arr[tbound][axis]; 
    int i = lbound - 1; 
    for(int j = lbound; j < tbound; j++){
        if(arr[j][axis] <= pivot){
	        i++;
            for(int q = 0; q < dimension + mode; q++){
	            skeleSwap(&arr[i][q], &arr[j][q]);
            }
	    }
    }
    for(int q = 0; q < dimension + mode; q++){
        skeleSwap(&arr[i + 1][q], &arr[tbound][q]);
    }
    int temp_pivot = i+1;
    return temp_pivot;
}

void skeleQuickSort(double **arr, int lbound, int tbound, int axis,int mode){
    //lbound => lower bounds
    //tbound => top bounds/length of array
    //Axis is axis to sort alng
    //1 => x
    //2 => y
    //3 => z
    if(lbound < tbound){
        double pivot = skelePartition(arr, lbound, tbound, axis, mode);
        //fprintf(stdout,"got pivot at%f\n",pivot);
        if(pivot - 1 - lbound > 1){
            skeleQuickSort(arr, lbound, pivot - 1, axis, mode);
        }
        if(tbound - pivot + 1 > 1){
            skeleQuickSort(arr, pivot + 1, tbound, axis, mode);
        }
    }
}


