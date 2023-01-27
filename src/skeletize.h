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
    //double *leftbound;
    //double *rightbound;
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
    int dim = *(&searchpoint + 1) - searchpoint;
    if(tree.flag == false){
        //Can go deeper, finds node to dive to 
        bool side = searchpoint[tree.axis] < tree.origin[tree.axis] ;
        if(side){
            retpoint = getNearest(searchpoint,*tree.left);
        }
        else{
            retpoint = getNearest(searchpoint,*tree.right);   
        }
        //next we will test if we need to go to the other side of the tree
        double axisnodedis = abs((int)(searchpoint[tree.axis] - tree.origin[tree.axis]));//This only measure along one axis
        if(axisnodedis < *retpoint[dim]){
            double **tempretpoint;
            if(side){
                tempretpoint = getNearest(searchpoint,*tree.right);
            }
            else{
                tempretpoint = getNearest(searchpoint,*tree.left);   
            }
            if(tempretpoint[dim] < retpoint[dim]){
                retpoint = tempretpoint;
            }
        }
        //finally check the node
        double nodedis = getDistance(searchpoint,tree.origin);
        if(nodedis < *retpoint[dim]){
            retpoint[0] = tree.origin[0];
            retpoint[1] = tree.origin[1];
            if(dim == 3){
                retpoint[2] = tree.origin[2];
                retpoint[3] = &nodedis;
            } 
            else{
                retpoint[2] = &nodedis;
            }
        }
    }
    else{
        //Lowest search level.
        double lowdis = 0.0;
        int length = *(&tree.origin + 1) - tree.origin;
        for(int i = 0; i < length; i++){
            double tdis = getDistance(searchpoint,&tree.origin[i]);
            if(lowdis == 0.0 || tdis < lowdis){
                lowdis = tdis;
                retpoint[0] = &tree.origin[i][0];
                retpoint[1] = &tree.origin[i][1];
                if(dim == 3){
                    retpoint[2] = &tree.origin[i][2];
                    retpoint[3] = &lowdis;
                } 
                else{
                    retpoint[2] = &lowdis;
                }
            }
        }
    }
    return retpoint;
}

void CreateStructure(double **points,struct kdleaf tree,int axis,int *length,int lower,int upper){//,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    //tree.leftbound = leftbound;
    //tree.rightbound = rightbound;
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
        CreateStructure(points,*tree.left,axis + 1,length,lower,nodeindex - 1);
        CreateStructure(points,*tree.right,axis+ 1,length,nodeindex + 1,upper);
    }
    else if(upper-lower > 0){
        tree.flag = true;
        *tree.origin = slice(*points,lower,upper);//will make origin now all 
    }
}
double getRadius(double **point,double **interface,int *dim){
    //point has norm data attached
    //dim is size of data(makes life easier)
    double top;
    double bot = getDistance(interface,point);
    for(int i = 0; i < *dim; i++){
        //goes though each dimension
        top += *point[*dim + i] * (*point[i] - *interface[i]);
    }
    double theta = acos(top / bot);
    double ret = (bot / (2 * cos(theta)));
    return ret;
}
double **makeSkeleton(double **points,struct kdleaf tree,int *length,int *max){
    //length is dim
    //max is size of list
    printf("starting process");
    double guessr = *max;
    double **skeleton;
    double **centerPoint;
    double *radius;
    double **interfacePoint;
    for(int i = 0; i < *max; i++){
        //Goes through each element of the list & generates a skeleton point
        //First step is to make an initial guess
        bool completeCase = true;
        int index = 0;
        double x = points[i][0] + guessr;
        double y = points[i][1] + guessr;
        double z;
        double **tpoint;
        tpoint[0] = &x;
        tpoint[1] = &y;
        if(*length == 3){
            double z = points[i][2] + guessr;
            tpoint[2] = &z;
        }
        //We get a interface point and calulate the raidus to begin
        interfacePoint[index] = *getNearest(tpoint,tree);
        radius[index] = getRadius(&points[i],&interfacePoint[index],length);
        while(completeCase){
            //first calculate new center point based on the last radius given
            if(*length == 2){
                x = points[i][0] - radius[index] * points[i][2];
                y = points[i][1] - radius[index] * points[i][3];
                tpoint[0] = &x;
                tpoint[1] = &y;
            }
            else{
                x = points[i][0] - radius[index] * points[i][3];
                y = points[i][1] - radius[index] * points[i][4];
                z = points[i][2] - radius[index] * points[i][5];
                tpoint[0] = &x;
                tpoint[1] = &y;
                tpoint[2] = &z;
            }
            centerPoint[index] = *tpoint;
            //then calculate the closest point
            interfacePoint[index + 1] = *getNearest(tpoint,tree);
            radius[index + 1] = getRadius(&points[i],&interfacePoint[index + 1],length);
            if(fabs(radius[index] - radius[index + 1]) < 0.001){
                //convergance conditions
                //our center point should remain the same
                for(int j = 0; j < *length; j++){
                    //Create skeleton point
                    skeleton[i][j] = centerPoint[index][j];
                }
                skeleton[i][*length] = radius[index];
                completeCase = false;
            }
            index += 1;
            if(index >= 50){
                break;
            }
        }
    }
    return skeleton;
}
void skeletize(double **points){
    //initial setup
    struct skeleData data;
    struct kdleaf tree;
    data.points = points;
    //sizes
    int length = (*(&points[0] + 1) - points[0]) / 2;//gets dimension of skeleton data, should give 2 or 3
    printf("length(dim) is: %d",length);
    int max = *(&points + 1) - points;//gets amount of points being put in
    printf("we have %d points",max);
    CreateStructure(points,tree,0,&length,0,max);//make kd-tree
    
    //Next we will skeletonize all the points in our list
    double **skeleton = makeSkeleton(points,tree,&length,&max);
    //print results
    for(int i = 0; i < max; i++){
        if(length == 2){
            printf("Skeleton %d: [%f,%f] Radius: %f",i,skeleton[i][0],skeleton[i][1],skeleton[i][2]);
        }
        else{
            printf("Skeleton %d: [%f,%f,%f] Radius: %f",i,skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]);
        }
    }
}
