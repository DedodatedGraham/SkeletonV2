#include <stdio.h>
#include <math.h>
#include <time.h>

void Swap(double *a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int Partition(double **arr, int lbound, int tbound, int axis, int *length){
    double pivot = arr[tbound][axis]; 
    int i = lbound - 1; 
    for(int j = lbound; j < tbound; j++){
        if(arr[j][axis] <= pivot){
	        i++;
            if(*length == 2){
	            Swap(&arr[i][axis    ], &arr[j][axis    ]);
	            Swap(&arr[i][axis - 1], &arr[j][axis - 1]);
	            Swap(&arr[i][axis - 2], &arr[j][axis - 2]);
	            Swap(&arr[i][axis - 3], &arr[j][axis - 3]);
            }
            else{
	            Swap(&arr[i][axis    ], &arr[j][axis    ]);
	            Swap(&arr[i][axis - 1], &arr[j][axis - 1]);
	            Swap(&arr[i][axis - 2], &arr[j][axis - 2]);
	            Swap(&arr[i][axis - 3], &arr[j][axis - 3]);
	            Swap(&arr[i][axis - 4], &arr[j][axis - 4]);
	            Swap(&arr[i][axis - 5], &arr[j][axis - 5]);
            }
	    }
    }
    if(*length == 2){
	    Swap(&arr[i + 1][axis    ], &arr[tbound][axis    ]);
	    Swap(&arr[i + 1][axis - 1], &arr[tbound][axis - 1]);
	    Swap(&arr[i + 1][axis - 2], &arr[tbound][axis - 2]);
	    Swap(&arr[i + 1][axis - 3], &arr[tbound][axis - 3]);
    }
    else{
	    Swap(&arr[i + 1][axis    ], &arr[tbound][axis    ]);
	    Swap(&arr[i + 1][axis - 1], &arr[tbound][axis - 1]);
	    Swap(&arr[i + 1][axis - 2], &arr[tbound][axis - 2]);
	    Swap(&arr[i + 1][axis - 3], &arr[tbound][axis - 3]);
	    Swap(&arr[i + 1][axis - 4], &arr[tbound][axis - 4]);
	    Swap(&arr[i + 1][axis - 5], &arr[tbound][axis - 5]);
    }
    int temp_pivot = i+1;
    return temp_pivot;
}

void QuickSort(double **arr, int lbound, int tbound, int axis, int *length){
    //lbound => lower bounds
    //tbound => top bounds
    //Axis is axis to sort alng
    //1 => x
    //2 => y
    //3 => z
    if(lbound < tbound){
        double pivot = Partition(arr, lbound, tbound, axis, length);
        QuickSort(arr, lbound, pivot - 1, axis, length);
        QuickSort(arr, pivot + 1, tbound, axis, length);
    }
}


