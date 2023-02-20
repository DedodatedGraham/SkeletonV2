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
    int axis;
    int leng;
    //if its the lowest layer, then it will contain a colletion of points
    //lower layers
    //double *leftbound;
    //double *rightbound;
    struct kdleaf *left,*right; 
};
//method for slicing
double **slice(double *array, int start, int end,int dim){
    fprintf(stdout,"slicing from %d to %d \n",start,end);
    int numElements = (end - start + 1);//find length of wanted slice
    int numBytes = sizeof(double) * (dim * 2) * numElements;//find byte size
    double **slice = (double**)malloc(numBytes);//allocate
    memcpy(slice, array + start, numBytes);//shift
    return slice;
}
//method for getting distance
double getDistance(double *point1,double *point2,int *dim){
    //sqrt sum of squares for distance
    double distance = 0;
    for(int i = 0; i < *dim; i++){
        double tsep = point1[i] - point2[i];
        double tabs = fabs(tsep);
        double tpow = pow(tabs,2);
        distance += tpow;
    }
    distance = sqrt(distance);
    return distance;
}

double *getNearest(double *searchpoint,struct kdleaf *kdstruct, int *length,int *dim,double *ignorepoint,double *lowestdistance){
    double *retpoint;
    if(kdstruct->flag){
        //Lowest search level.
        double lowdis = 0.0;
        for(int i = 0; i < kdstruct->leng; i++){
            double tdis = getDistance(searchpoint,kdstruct->origin[i],dim);
            if((lowdis == 0.0 || tdis < lowdis) && getDistance(ignorepoint,kdstruct->origin[i],dim) != 0.0){
                lowdis = tdis;
                retpoint = kdstruct->origin[i];
                *lowestdistance = lowdis;
            }
        }
    }
    else{
        //can go deeper
        bool side = searchpoint[kdstruct->axis] < kdstruct->origin[0][kdstruct->axis] ;
        if(side){
            retpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,lowestdistance);
        }
        else{
            retpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,lowestdistance);   
        }
        //next we will test if we need to go to the other side of the struct
        double axisnodedis = fabs((int)(searchpoint[kdstruct->axis] - kdstruct->origin[0][kdstruct->axis]));//This only measure along one axis
        if(axisnodedis < retpoint[*dim]){
            double *tempretpoint;
            if(side){
                tempretpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,lowestdistance);
            }
            else{
                tempretpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,lowestdistance);   
            }
            if(tempretpoint[*dim] < retpoint[*dim]){
                retpoint = tempretpoint;
            }
        }
        //finally check the node
        double nodedis = getDistance(searchpoint,kdstruct->origin[0],dim);
        if(getDistance(ignorepoint,kdstruct->origin[0],dim) != 0.0 && nodedis < *lowestdistance){
            retpoint = kdstruct->origin[0];
            *lowestdistance = nodedis;
        }
    }
    return retpoint;
}


//Destroy
void kdDestroy(struct kdleaf *kdstruct){
    free(kdstruct->origin);
    kdstruct->origin = NULL;
    if(!kdstruct->flag){
        kdDestroy(kdstruct->left);
        kdDestroy(kdstruct->right);
    }
    free(kdstruct);
    kdstruct->origin = NULL;
    kdstruct->flag = NULL;
    kdstruct = NULL;
}
struct kdleaf *createLeaf(int axis, double **points){
    struct kdleaf *kdstruct = (struct kdleaf*)malloc(sizeof(struct kdleaf));
    kdstruct->axis = axis;
    kdstruct->origin = points;
    return kdstruct;
}
struct kdleaf *CreateStructure(double **points,struct kdleaf **headkdstruct,int axis,int *length,int lower,int upper){//,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    if(axis == *length){
        axis = 0;//will correct back to right axis
    }
    struct kdleaf *currentkdstruct = *headkdstruct;
    if(upper-lower > 5){
        skeleQuickSort(points,lower,upper,axis,length);
        int nodeindex = (upper + lower) / 2; 
        double **nodepoint = (double**)malloc(sizeof(double*));
        nodepoint[0] = points[nodeindex];
        currentkdstruct = createLeaf(axis,nodepoint);
        bool flag = false;
        currentkdstruct->flag = flag; 
        
        //points now sorted, Next applying our node point and split the list to the appropiate list

        currentkdstruct->left = CreateStructure(points,&currentkdstruct->left,axis + 1,length,lower,nodeindex - 1);
        currentkdstruct->right = CreateStructure(points,&currentkdstruct->right,axis+ 1,length,nodeindex + 1,upper);
        *headkdstruct = currentkdstruct;
    }
    else if(upper-lower > 0){
        double **hold = (double**)malloc((upper - lower) * sizeof(double*));
        for(int i = 0; i < upper - lower; i++){
            hold[i] = points[lower + i];
        }
        currentkdstruct = createLeaf(axis,hold);
        bool flag = true;
        currentkdstruct->flag = flag;
        int leng = upper-lower;
        currentkdstruct->leng = leng;
    }
    return currentkdstruct;
}
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
double getRadius(double *point,double *interface,int *dim){
    //point has norm data attached
    //dim is size of data(makes life easier)
    double top1 = 0.0;
    double bot = getDistance(interface,point,dim);
    for(int i = 0; i < *dim; i++){
        //goes though each dimension
        double diff1 = point[i] - interface[i];
        double add1 = point[*dim+i] * diff1;
        top1 += add1;
    }
    double inside1 = top1/bot; 
    if(inside1 > 1){
        inside1 = 1;
    }
    if(inside1 < -1){
        inside1 = -1;
    }
    double theta1 = acos(inside1);
    double ret1 = (bot / (2 * cos(theta1)));
    if(ret1 > 0){
        return ret1;
    }
    else{
        fprintf(stdout,"reversed radius\n");
        return -ret1;
    }
}
void outputskeleton(double *points, int *dim, char path[80]){
    //we save each point as it comes in
    FILE *fp = fopen(path,"a");
    if(*dim == 2){
        fprintf(fp,"%f %f %f\n",points[0],points[1],points[2]);
    }
    else{
        fprintf(fp,"%f %f %f %f\n",points[0],points[1],points[2],points[3]);
    }
    fflush(fp);
    fclose(fp);
}
void makeSkeleton(double **points,struct kdleaf *kdstruct,int *dim,int *length,double *mindis,char path[80]){
    int MAXCYCLES = 50;

    //fprintf(stdout,"starting process\n");
    double guessr = *length;
    double **skeleton = (double**)malloc((*length + 1) * sizeof(double*));
    double **centerPoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    double *radius = (double*)malloc(MAXCYCLES * sizeof(double));
    double **interfacePoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    for(int i = 0; i < *length + 1; i++){
        //fprintf(stdout,"\ncalculating point: %d\n",i);
        //Goes through each element of the list & generates a skeleton point
        //First step is to make an initial guess
        bool completeCase = true;
        int index = 0;
        double *ttpoint = (double*)malloc(*dim * sizeof(double));
        double *ignorepoint = (double*)malloc(*dim * sizeof(double));
        if(*dim == 2){
            double x = points[i][0] - points[i][2] * guessr;
            double y = points[i][1] - points[i][3] * guessr;
            double ix = points[i][0];
            double iy = points[i][1];
            ttpoint[0] = x;
            ttpoint[1] = y;
            ignorepoint[0] = ix;
            ignorepoint[1] = iy;
            //ttpoint = NULL;
            //ignorepoint = NULL;
            //fprintf(stdout,"\n(%d,%d)\n",i,index);
            //fprintf(stdout,"x = %f or %f\n",ipoint[0],points[i][0]);
            //fprintf(stdout,"y = %f or %f\n",ipoint[1],points[i][1]);
        }
        else{
            double x = points[i][0] - points[i][3] * guessr;
            double y = points[i][1] - points[i][4] * guessr;
            double z = points[i][2] - points[i][5] * guessr;
            double ix = points[i][0];
            double iy = points[i][1];
            double iz = points[i][2];
            double *ttpoint = (double*)malloc(*dim * sizeof(double));
            double *ignorepoint = (double*)malloc(*dim * sizeof(double));
            ttpoint[0] = x;
            ttpoint[1] = y;
            ttpoint[2] = z;
            ignorepoint[0] = ix;
            ignorepoint[1] = iy;
            ignorepoint[2] = iz;
            //free(ttpoint);
            //free(ignorepoint);
            //ttpoint = NULL;
            //ignorepoint = NULL;
        }
        double lowestdistance = 0; 
        centerPoint[index] = ttpoint;
        interfacePoint[index] = getNearest(ttpoint,kdstruct,length,dim,ignorepoint,&lowestdistance);
        fprintf(stdout,"\nskeleton for point%d x:%f y:%f nx:%f ny:%f\n",i,points[i][0],points[i][1],points[i][2],points[i][3]);
        fprintf(stdout,"center%d x:%f y:%f \n",index,ttpoint[0],ttpoint[1]);
        fprintf(stdout,"interface point%d x:%f y:%f\n",index,interfacePoint[index][0],interfacePoint[index][1]); 
        double* sendpoint = points[i];
        radius[index] = getRadius(sendpoint,interfacePoint[index],dim);
        fprintf(stdout,"rad (%d,%d) = %f\n",i,index,radius[index]);
        while(completeCase){
            fprintf(stdout,"calculating:%d\n",index+1);
            //first calculate new center point based on the last radius given
            //fprintf(stdout,"dim = %d\n",*dim);
            if(*dim == 2){
                double x = points[i][0] - points[i][2] * radius[index];
                double y = points[i][1] - points[i][3] * radius[index];
                double ix = points[i][0];
                double iy = points[i][1];
                ttpoint[0] = x;
                ttpoint[1] = y;
                ignorepoint[0] = ix;
                ignorepoint[1] = iy;
                //free(ttpoint);
                //free(ignorepoint);
                //ttpoint = NULL;
                //ignorepoint = NULL;
                //fprintf(stdout,"(%d,%d)\n",i,index+1);
                //fprintf(stdout,"x = %f or %f\n",ipoint[0],points[i][0]);
                //fprintf(stdout,"y = %f or %f\n",ipoint[1],points[i][1]);
            }
            else{
                double x = points[i][0] - points[i][3] * radius[index];
                double y = points[i][1] - points[i][4] * radius[index];
                double z = points[i][2] - points[i][5] * radius[index];
                double ix = points[i][0];
                double iy = points[i][1];
                double iz = points[i][2];
                double *ttpoint = (double*)malloc(*dim * sizeof(double));
                double *ignorepoint = (double*)malloc(*dim * sizeof(double));
                ttpoint[0] = x;
                ttpoint[1] = y;
                ttpoint[2] = z;
                ignorepoint[0] = ix;
                ignorepoint[1] = iy;
                ignorepoint[2] = iz;
                //tpoint = ttpoint;
                //ipoint = ignorepoint;
                //free(ttpoint);
                //free(ignorepoint);
                //ttpoint = NULL;
                //ignorepoint = NULL;
            }
            centerPoint[index] = ttpoint;
            fprintf(stdout,"center%d x:%f y:%f \n",index + 1,ttpoint[0],ttpoint[1]);
            interfacePoint[index + 1] = getNearest(centerPoint[index],kdstruct,length,dim,ignorepoint,&lowestdistance);
            fprintf(stdout,"interface point%d x:%f y:%f\n",index + 1,interfacePoint[index + 1][0],interfacePoint[index + 1][1]); 
            radius[index + 1] = getRadius(points[i],interfacePoint[index + 1],dim);
            fprintf(stdout,"rad (%d,%d) = %f\n",i,index+1,radius[index+1]);
            double distancecomp = getDistance(interfacePoint[index + 1],points[i],dim);
            if(radius[index] != 0. && fabs(radius[index] - radius[index + 1]) < *mindis){
                //convergance conditions
                //our center point should remain the same
                skeleton[i] = centerPoint[index];
                skeleton[i][*dim] = radius[index + 1];
                fprintf(stdout,"converge skelept:[x=%f,y=%f,r=%f]\n",skeleton[i][0],skeleton[i][1],skeleton[i][2]);
                fprintf(stdout,"point : [x=%f,y=%f]\n",points[i][0],points[i][1]);
                fprintf(stdout,"interfacepoint : [x=%f,y=%f]\n",interfacePoint[index + 1][0],interfacePoint[index + 1][1]);
                outputskeleton(skeleton[i],dim,path);
                completeCase = false;
            }
            else if(index > 0 && distancecomp < radius[index + 1]){
                skeleton[i] = centerPoint[index-1];
                skeleton[i][*dim] = radius[index];
                //fprintf(stdout,"\n(%d,%d)\n",i,index);
                fprintf(stdout,"discomp skelept:[x=%f,y=%f,r=%f]\n",skeleton[i][0],skeleton[i][1],skeleton[i][2]);
                fprintf(stdout,"point : [x=%f,y=%f]\n",points[i][0],points[i][1]);
                //fprintf(stdout,"centpoint-2 : [x=%f,y=%f],%f\n",centerPoint[index - 2][0],centerPoint[index - 2][1],radius[index-1]);
                //fprintf(stdout,"centpoint-1 : [x=%f,y=%f],%f\n",centerPoint[index - 1][0],centerPoint[index - 1][1],radius[index]);
                //fprintf(stdout,"centpoint-0 : [x=%f,y=%f],%f\n",centerPoint[index][0],centerPoint[index][1],radius[index+1]);
                fprintf(stdout,"interfacepoint : [x=%f,y=%f]\n",interfacePoint[index + 1][0],interfacePoint[index + 1][1]);
                outputskeleton(skeleton[i],dim,path);
                completeCase = false;
            }
            //else if(distancecomp < *mindis * 10){
            //    skeleton[i] = centerPoint[index - 1];
            //    skeleton[i][*dim] = radius[index];
            //    fprintf(stdout,"skelept:[x=%f,y=%f,r=%f]\n",skeleton[i][0],skeleton[i][1],skeleton[i][2]);
            //    completeCase = false;
            //}
            index += 1;
            if(index >= MAXCYCLES){
                fprintf(stdout,"broken\n");
                break;
            }
        }
        free(ttpoint);
        free(ignorepoint);
    }
    free(radius);
    free(centerPoint);
    free(interfacePoint);
    free(skeleton);
    radius = NULL;
    centerPoint = NULL;
    interfacePoint = NULL;
    skeleton = NULL;
}
void skeletize(double **points,int *length,int *dim,char path[80],double *mindis){
    //initial setup
    struct kdleaf *kdstruct;
    CreateStructure(points,&kdstruct,0,dim,0,*length);//make kd-struct
    makeSkeleton(points,kdstruct,dim,length,mindis,path);
    for(int i = 0;i < *length + 1; i++){
        free(points[i]);
    }
    free(points);
    points = NULL;
    kdDestroy(kdstruct);
}
 
