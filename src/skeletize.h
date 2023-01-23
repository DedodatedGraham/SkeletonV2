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
    int axis;
    //if its the lowest layer, then it will contain a colletion of points
    //lower layers
    double *leftbound;
    double *rightbound;
    struct kdleaf *left,*right; 
};
//method for slicing
double *slice(double *array, int start, int end){
    int numElements = (end - start + 1);//find length of wanted slice
    int numBytes = sizeof(int) * numElements;//find byte size

    double *slice = malloc(numBytes);//allocate
    memcpy(slice, array + start, numBytes);//shift

    return slice;
}
//method for getting distance
double getDistance(double **point1,double **point2){
    //sqrt sum of squares for distance
    int length = *(&point1 + 1) - point1;
    double distance = 0;
    for(int i = 0; i < length; i++){
        distance += pow(abs((int)(point1[i] - point2[i])),2);
    }
    distance = sqrt(distance);
    return distance;
}
double **getNearest(double **searchpoint,struct kdleaf tree){
    double ** retpoint;
    if(tree.flag == false){
        //Can go deeper, finds node to dive to 
        if(searchpoint[tree.axis] < tree.origin[tree.axis]){
            retpoint = getNearest(searchpoint,*tree.left);
        }
        else{
            retpoint = getNearest(searchpoint,*tree.right);   
        }
        //next we will test distances and bounds and see if we should just return

    }
    else{
        //Lowest search level.
        double lowdis = 0.0;
        int length = *(&tree.origin + 1) - tree.origin;
        for(int i = 0; i < length; i++){
            double tdis = getDistance(searchpoint,&tree.origin[i]);
            if(lowdis == 0.0 || tdis < lowdis){
                lowdis = tdis;
                retpoint = &tree.origin[i];
            }
        }
    }
    return retpoint;
}

void CreateStructure(double **points,struct kdleaf tree,int axis,int *length,int lower,int upper,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    tree.leftbound = leftbound;
    tree.rightbound = rightbound;
    tree.axis = axis;
    if(upper-lower > 5){
        tree.flag = false;
        QuickSort(points,lower,upper,axis,length);
        //points now sorted, Next applying our node point and split the list to the appropiate list
        int nodeindex = (upper - lower) / 2; 
        tree.origin = &points[nodeindex];
        //getting new bounds for lower layers
        //bounds change depending on axis dividing
        
        if(axis == *length - 1){
            axis = -1;//will correct back to right axis
        }
        
        if(*length == 2){
            double *newrightbound;
            double *newleftbound;
            if(axis == 0){
                //x axis
                newrightbound[0] = *tree.origin[0];
                newrightbound[1] = rightbound[1];
                newleftbound[0] = *tree.origin[0];
                newleftbound[1] = leftbound[1];
            }
            else{
                //yaxis
                newrightbound[0] = rightbound[0];
                newrightbound[1] = *tree.origin[1];
                newleftbound[0] = leftbound[0];
                newleftbound[1] = *tree.origin[1];
            }
            CreateStructure(points, *tree.left, axis + 1, length, lower, nodeindex - 1,leftbound,newrightbound);
            CreateStructure(points, *tree.right,axis + 1, length, nodeindex + 1, upper,newleftbound,rightbound);
        }
        else{
            double *newrightbound;
            double *newleftbound;
            if(axis == 0){
                //x axis
                newrightbound[0] = *tree.origin[0];
                newrightbound[1] = rightbound[1];
                newrightbound[2] = rightbound[2];
                newleftbound[0] = *tree.origin[0];
                newleftbound[1] = leftbound[1];
                newleftbound[2] = leftbound[2];
            }
            else if(axis == 1){
                //y axis
                newrightbound[0] = rightbound[0];
                newrightbound[1] = *tree.origin[1];
                newrightbound[2] = rightbound[2];
                newleftbound[0] = leftbound[0];
                newleftbound[1] = *tree.origin[1];
                newleftbound[2] = leftbound[2];
            }
            else{
                //z axis
                newrightbound[0] = rightbound[0];
                newrightbound[1] = rightbound[1];
                newrightbound[2] = *tree.origin[2];
                newleftbound[0] = leftbound[0];
                newleftbound[1] = leftbound[1];
                newleftbound[2] = *tree.origin[2];
            }
            CreateStructure(points, *tree.left, axis + 1, length, lower, nodeindex - 1,leftbound,newrightbound);
            CreateStructure(points, *tree.right,axis + 1, length, nodeindex + 1, upper,newleftbound,rightbound);
        }
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
    double *leftbound;
    double *rightbound;
    if(length == 2){
        leftbound[0] = xbounds[0];
        leftbound[1] = ybounds[0];
        rightbound[0] = xbounds[1];
        rightbound[1] = ybounds[1];
    }
    else{
        leftbound[0] = xbounds[0];
        leftbound[1] = ybounds[0];
        leftbound[2] = zbounds[0];
        rightbound[0] = xbounds[1];
        rightbound[1] = ybounds[1];
        rightbound[2] = zbounds[1];
    }
    //Builds tree
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
