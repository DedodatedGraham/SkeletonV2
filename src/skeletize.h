#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "skeleQuicksort.h"

//First we define a data structure, it will point to the data points
struct skeleData{
    double **points;//This is a pointer to the list of points and normals
    //If 2D then [x, y, normx, normy]
    //if 3D then [x, y, z, normx, normy, normz]
    double threshold;//This is a a double which is a distance average of the points
};

struct kdleaf{
    //Bounds and position
    double **origin;//node point [x,y] or [x,y,z]
    bool flag;
    int depth;
    //if its the lowest layer, then it will contain a colletion of points
    //lower layers
    double *leftbound;
    double *rightbound;
    struct kdleaf *left,*right; 
};

double *slice(double *array, int start, int end) {
    int numElements = (end - start + 1);//find length of wanted slice
    int numBytes = sizeof(int) * numElements;//find byte size

    double *slice = malloc(numBytes);//allocate
    memcpy(slice, array + start, numBytes);//shift

    return slice;
}
double **getNearest(double *searchpoint,struct kdleaf tree,int axis){
    if(tree.flag == false){
        //Can go deeper
    }
    else{
        //Lowest search level
    }
}

void CreateStructure(double **points,struct kdleaf tree,int axis,int *length,int lower,int upper,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    tree.leftbound = leftbound;
    tree.rightbound = rightbound;
    if(upper-lower > 5){
        tree.flag = false;
        QuickSort(points,lower,upper,axis,length);
        //points now sorted, Next applying our node point and split the list to the appropiate list
        int nodeindex = (upper - lower) / 2; 
        tree.origin = &points[nodeindex];
        //getting new bounds for lower layers
        //bounds change depending on axis dividing
        if(axis == 2){
            //3D default
        } 
        else if(axis == 1){
            //division in Y
            
        }
        else{
            //division in x
            if(*length == 2){
                double *newrightbound [2] = {tree.origin[0],&rightbound[1]};
            }
            else{

            }
        }
        

        //then we spilt the arrays and go deeper
        if(axis == *length - 1){
            axis = -1;//will correct back to right axis
        }
        CreateStructure(points, *tree.left, axis + 1, length, lower, nodeindex - 1,leftbound,newrightbound);
        CreateStructure(points, *tree.right,axis + 1, length, nodeindex + 1, upper,newleftbound,rightbound);
    }
    else if(upper-lower > 0){
        tree.flag = true;
        *tree.origin = slice(*points,lower,upper);//will make origin now all 
    }
}

void skeltize(double **points){
    //initial setup
    struct skeleData data;
    struct kdleaf tree;
    data.points = points;
    //sizes
    int length = (*(&points[0] + 1) - points[0]) / 2;//gets dimension of skeleton data, should give 2 or 3
    printf("length(dim) is: %d",length);
    int max = *(&points + 1) - points;//gets amount of points being put in
    printf("we have %d points",max);
    //gets bounds
    double xbounds [2];
    double ybounds [2];
    double zbounds [2];
    xbounds[0] = points[0][0];
    xbounds[1] = points[1][0];
    ybounds[0] = points[0][1];
    ybounds[1] = points[1][1];
    if(length == 3){
        zbounds[0] = points[0][2];
        zbounds[1] = points[1][2];
    }
    for(int i = 2; i < max; i ++){ 
        if(points[i][0] < xbounds[0]){
            xbounds[0] = points[i][0];
        }
        else if(points[i][0] > xbounds[1]){
            xbounds[1] = points[i][0];
        }
        if(points[i][1] < ybounds[0]){
            ybounds[0] = points[i][1];
        }
        else if(points[i][1] > ybounds[1]){
            ybounds[1] = points[i][1];
        }
        if(length == 3){
            if(points[i][2] < zbounds[0]){
                zbounds[0] = points[i][2];
            }
            else if(points[i][2] > zbounds[1]){
                zbounds[1] = points[i][2];
            }
        }
    }
    if(length == 2){
        double *leftbound [2] = {&xbounds[0],&ybounds[0]};
        double *rightbound [2] = {&xbounds[1],&ybounds[1]};
    }
    else{
        double *leftbound [3] = {&xbounds[0],&ybounds[0],&zbounds[0]};
        double *rightbound [3] = {&xbounds[1],&ybounds[1],&zbounds[0]};
    }
    CreateStructure(points,tree,0,&length,0,max,leftbound,rightbound);//make kd-tree

}




//void Swap(double *a, double *b){
//    double temp = *a;
//    *a = *b;
//    *b = temp;
//}
//
//int Partition(double **arr, int lbound, int tbound, int index){
//    //lbound => lower bounds
//    //tbound => top bounds
//    double pivot = arr[tbound][index]; 
//    int i = lbound - 1; 
//
//    for(int j = lbound; j < tbound; j++){
//        if(arr[j][index] <= pivot){
//	    i++;
//	    Swap(&arr[i][index], &arr[j][index]);
//	    Swap(&arr[i][index - 1], &arr[j][index - 1]);
//	    Swap(&arr[i][index - 2], &arr[j][index - 2]);
//	    }
//    }
//    Swap(&arr[i + 1][index], &arr[tbound][index]);
//    Swap(&arr[i + 1][index - 1], &arr[tbound][index - 1]);
//    Swap(&arr[i + 1][index - 2], &arr[tbound][index - 2]);
//    int temp_pivot = i+1;
//return temp_pivot;
//}
//
//void QuickSort(double **arr, int lbound, int tbound, int index){
//    if(lbound < tbound){
//        double pivot = Partition(arr, lbound, tbound, index);
//        QuickSort(arr, lbound, pivot - 1, index);
//        QuickSort(arr, pivot + 1, tbound, index);
//    }
//}
