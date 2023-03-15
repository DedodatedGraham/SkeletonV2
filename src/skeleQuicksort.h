#include <stdio.h>
#include <math.h>
#include <time.h>

void skeleSwap(double *a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int skelePartition(double **arr, int lbound, int tbound, int axis, int *length){
    double pivot = arr[tbound][axis]; 
    int i = lbound - 1; 
    for(int j = lbound; j < tbound; j++){
        if(arr[j][axis] <= pivot){
	        i++;
            if(*length == 2){
	            skeleSwap(&arr[i][0], &arr[j][0]);
	            skeleSwap(&arr[i][1], &arr[j][1]);
	            skeleSwap(&arr[i][2], &arr[j][2]);
	            skeleSwap(&arr[i][3], &arr[j][3]);
            }
            else{
	            skeleSwap(&arr[i][0], &arr[j][0]);
	            skeleSwap(&arr[i][1], &arr[j][1]);
	            skeleSwap(&arr[i][2], &arr[j][2]);
	            skeleSwap(&arr[i][3], &arr[j][3]);
	            skeleSwap(&arr[i][4], &arr[j][4]);
	            skeleSwap(&arr[i][5], &arr[j][5]);
            }
	    }
    }
    if(*length == 2){
	    skeleSwap(&arr[i + 1][0], &arr[tbound][0]);
	    skeleSwap(&arr[i + 1][1], &arr[tbound][1]);
	    skeleSwap(&arr[i + 1][2], &arr[tbound][2]);
	    skeleSwap(&arr[i + 1][3], &arr[tbound][3]);
    }
    else{
	    skeleSwap(&arr[i + 1][0], &arr[tbound][0]);
	    skeleSwap(&arr[i + 1][1], &arr[tbound][1]);
	    skeleSwap(&arr[i + 1][2], &arr[tbound][2]);
	    skeleSwap(&arr[i + 1][3], &arr[tbound][3]);
	    skeleSwap(&arr[i + 1][4], &arr[tbound][4]);
	    skeleSwap(&arr[i + 1][5], &arr[tbound][5]);
    }
    int temp_pivot = i+1;
    return temp_pivot;
}

void skeleQuickSort(double **arr, int lbound, int tbound, int axis, int *length){
    //lbound => lower bounds
    //tbound => top bounds
    //Axis is axis to sort alng
    //1 => x
    //2 => y
    //3 => z
    if(lbound < tbound){
        double pivot = skelePartition(arr, lbound, tbound, axis, length);
        //fprintf(stdout,"got pivot at%f\n",pivot);
        if(pivot - 1 - lbound > 1){
            skeleQuickSort(arr, lbound, pivot - 1, axis, length);
        }
        if(tbound - pivot + 1 > 1){
            skeleQuickSort(arr, pivot + 1, tbound, axis, length);
        }
    }
}


