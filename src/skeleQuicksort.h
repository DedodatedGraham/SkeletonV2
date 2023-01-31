#include <stdio.h>
#include <math.h>
#include <time.h>

void skeleSwap(double *a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int skelePartition(double **arr, int lbound, int tbound, int axis, int *length){
    //fprintf(stdout,"partitioning %d, on axis %d\n",tbound,axis);
    double pivot = arr[tbound][axis]; 
    //fprintf(stdout,"piv = %f\n",pivot);
    int i = lbound - 1; 
    for(int j = lbound; j < tbound; j++){
        if(arr[j][axis] <= pivot){
	        i++;
            //fprintf(stdout,"swapping %d & %d\n",i,j);
            if(*length == 2){
	            skeleSwap(&arr[i][axis    ], &arr[j][axis    ]);
	            skeleSwap(&arr[i][axis - 1], &arr[j][axis - 1]);
	            skeleSwap(&arr[i][axis - 2], &arr[j][axis - 2]);
	            skeleSwap(&arr[i][axis - 3], &arr[j][axis - 3]);
            }
            else{
	            skeleSwap(&arr[i][axis    ], &arr[j][axis    ]);
	            skeleSwap(&arr[i][axis - 1], &arr[j][axis - 1]);
	            skeleSwap(&arr[i][axis - 2], &arr[j][axis - 2]);
	            skeleSwap(&arr[i][axis - 3], &arr[j][axis - 3]);
	            skeleSwap(&arr[i][axis - 4], &arr[j][axis - 4]);
	            skeleSwap(&arr[i][axis - 5], &arr[j][axis - 5]);
            }
	    }
    }
    //fprintf(stdout,"swapping %d & %d\n",i + 1,tbound);
    if(*length == 2){
	    skeleSwap(&arr[i + 1][axis    ], &arr[tbound][axis    ]);
	    skeleSwap(&arr[i + 1][axis - 1], &arr[tbound][axis - 1]);
	    skeleSwap(&arr[i + 1][axis - 2], &arr[tbound][axis - 2]);
	    skeleSwap(&arr[i + 1][axis - 3], &arr[tbound][axis - 3]);
    }
    else{
	    skeleSwap(&arr[i + 1][axis    ], &arr[tbound][axis    ]);
	    skeleSwap(&arr[i + 1][axis - 1], &arr[tbound][axis - 1]);
	    skeleSwap(&arr[i + 1][axis - 2], &arr[tbound][axis - 2]);
	    skeleSwap(&arr[i + 1][axis - 3], &arr[tbound][axis - 3]);
	    skeleSwap(&arr[i + 1][axis - 4], &arr[tbound][axis - 4]);
	    skeleSwap(&arr[i + 1][axis - 5], &arr[tbound][axis - 5]);
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
    //fprintf(stdout,"starting sort %d %d\n",lbound,tbound);
    if(lbound < tbound){
        double pivot = skelePartition(arr, lbound, tbound, axis, length);
        //fprintf(stdout,"got pivot at%f\n",pivot);
        if(pivot - 1 - lbound > 1){
            skeleQuickSort(arr, lbound, pivot - 1, axis, length);
        }
        if(tbound - pivot + 1 > 1){
            skeleQuickSort(arr, pivot + 1, tbound, axis, length);
        }
        //fprintf(stdout,"sort completed\n");
    }
    //else{
    //    fprintf(stdout,"sort skipped\n");
    //}
}


