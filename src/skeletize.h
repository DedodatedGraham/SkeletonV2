#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "skeleQuicksort.h"

//First we define a data structure, it will point to the data points
struct kdleaf{
    //Bounds and position
    double **origin;//node point [x,y] or [x,y,z] @ pos 0
                    //If lowest layer, then list of points
    bool flag;
    int axis;
    int leng;
    struct kdleaf *left,*right; 
};

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

double *getNearest(double *searchpoint,struct kdleaf *kdstruct, int *length,int *dim,double **ignorepoint,int *ileng,double *lowestdistance){
    double *retpoint;
    if(kdstruct->flag){
        //fprintf(stdout,"lowest level\n");
        //Lowest search level.
        double lowdis = 0.0;
        for(int i = 0; i < kdstruct->leng; i++){
            //fprintf(stdout,"\n");
            //fprintf(stdout,"leng kdpoints = %d\n",kdstruct->leng);
            double tdis = getDistance(searchpoint,kdstruct->origin[i],dim);
            if(*ileng == 1){
                if((lowdis == 0.0 || tdis < lowdis) && getDistance(ignorepoint[0],kdstruct->origin[i],dim) != 0.0){
                    lowdis = tdis;
                    retpoint = kdstruct->origin[i];
                    *lowestdistance = lowdis;
                }
            }
            else{
                bool passes = true;
                for(int j = 0; j < *ileng; j++){
                    //fprintf(stdout,"\n(%d,%d)\n",j,*ileng - 1);
                    //fprintf(stdout,"[%f,%f]\n",kdstruct->origin[i][0],kdstruct->origin[i][1]);
                    //fprintf(stdout,"[%f,%f]\n",ignorepoint[j][0],ignorepoint[j][1]);
                    if(kdstruct->origin[i] == ignorepoint[j]){
                        //fprintf(stdout,"broken\n");
                        passes = false;
                        break;
                    }
                }
                if((lowdis == 0.0 || tdis < lowdis) && passes){
                    //fprintf(stdout,"(%d/%d)passed\n",i,kdstruct->leng - 1);
                    lowdis = tdis;
                    retpoint = kdstruct->origin[i];
                    *lowestdistance = lowdis;
                    //fprintf(stdout,"distance = %f from point [%f,%f]\n",*lowestdistance,retpoint[0],retpoint[1]);
                }
                //else{
                //    //fprintf(stdout,"(%d/%d)didnt pass\n",i,kdstruct->leng - 1);
                //}
            }
            //fprintf(stdout,"got:%f\n",*lowestdistance);
        }
        //if(retpoint != NULL){
        //    //Here we double check the retpoint isnt bad
        //    if(getDistance(searchpoint,retpoint,dim) != *lowestdistance){
        //        fprintf(stdout,"bad point!\n");
        //        retpoint = NULL;
        //        *lowestdistance = 0.0;
        //    }
        //}
    }
    else{
        //can go deeper
        bool side = searchpoint[kdstruct->axis] < kdstruct->origin[0][kdstruct->axis];
        if(side){
            retpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,ileng,lowestdistance);
        }
        else{
            retpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,ileng,lowestdistance);   
        }
        //next we will test if we need to go to the other side of the struct
        double axisnodedis = fabs(searchpoint[kdstruct->axis] - kdstruct->origin[0][kdstruct->axis]);//This only measure along one axis
        if(axisnodedis < *lowestdistance || *lowestdistance == 0.0){
            double *tempretpoint;
            //double *templowdis = (double*)malloc(sizeof(double));
            double templowdis = 0.0;
            if(side){
                tempretpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,ileng,&templowdis);
            }
            else{
                tempretpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,ileng,&templowdis);   
            }
            //fprintf(stdout,"temp dist %f\n",templowdis);
            //fprintf(stdout,"temp point -> [%f,%f]\n",tempretpoint[0],tempretpoint[1]);
            if((*lowestdistance > templowdis|| *lowestdistance == 0.0) && templowdis != 0.0){
                retpoint = tempretpoint;
                *lowestdistance = templowdis;
            }
            //free(templowdis);
        }
        //finally check the node
        double nodedis = getDistance(searchpoint,kdstruct->origin[0],dim);
        if(*ileng == 1){
            if(getDistance(ignorepoint[0],kdstruct->origin[0],dim) != 0.0 && (nodedis < *lowestdistance || *lowestdistance == 0.0)){
                retpoint = kdstruct->origin[0];
                *lowestdistance = nodedis;
            }
        }
        else{
            bool passes = true;
            for(int j = 0; j < *ileng; j++){
                if(kdstruct->origin[0] == ignorepoint[j]){
                    passes = false;
                    break;
                }
            }
            if((nodedis < *lowestdistance || *lowestdistance == 0.0) && passes){
                retpoint = kdstruct->origin[0];
                *lowestdistance = nodedis;
                //fprintf(stdout,"set distance:%f\n",*lowestdistance);
            }
        }
    }
    //if(*lowestdistance != 0){
    //    fprintf(stdout,"lowest distance:%f\n",*lowestdistance);
    //}
    return retpoint;
}

////Here we have the Regional serach implemented for a kdtree
////Input the tree, output points inside - ignore points(multiple from stack if wanted)
//double** findBox(struct kdleaf *kdstruct,double *length,double *dim,double **box,double **ignorePoints){
//    double **returnPts; 
//    if(kdstruct->flag){
//        //Lowest search level.
//        //double lowdis = 0.0;
//        for(int i = 0; i < kdstruct->leng; i++){
//            //double tdis = getDistance(searchpoint,kdstruct->origin[i],dim);
//            //if((lowdis == 0.0 || tdis < lowdis) && getDistance(ignorepoint,kdstruct->origin[i],dim) != 0.0){
//            //    lowdis = tdis;
//            //    retpoint = kdstruct->origin[i];
//            //    *lowestdistance = lowdis;
//            //}
//            
//        }
//    }
//    else{
//        //can go deeper
//        bool side = searchpoint[kdstruct->axis] < kdstruct->origin[0][kdstruct->axis] ;
//        if(side){
//            retpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,lowestdistance);
//        }
//        else{
//            retpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,lowestdistance);   
//        }
//        //next we will test if we need to go to the other side of the struct
//        double axisnodedis = fabs(searchpoint[kdstruct->axis] - kdstruct->origin[0][kdstruct->axis]);//This only measure along one axis
//        if(axisnodedis < *lowestdistance){
//            double *tempretpoint;
//            double *templowdis = (double*)malloc(sizeof(double));
//            if(side){
//                tempretpoint = getNearest(searchpoint,kdstruct->right,length,dim,ignorepoint,templowdis);
//            }
//            else{
//                tempretpoint = getNearest(searchpoint,kdstruct->left,length,dim,ignorepoint,templowdis);   
//            }
//            if(*lowestdistance > *templowdis){
//                retpoint = tempretpoint;
//                *lowestdistance = *templowdis;
//            }
//            free(templowdis);
//        }
//        //finally check the node
//        double nodedis = getDistance(searchpoint,kdstruct->origin[0],dim);
//        if(getDistance(ignorepoint,kdstruct->origin[0],dim) != 0.0 && nodedis < *lowestdistance){
//            retpoint = kdstruct->origin[0];
//            *lowestdistance = nodedis;
//        }
//    }
//
//}



//Destroy tree
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
//make leaf of kd-tree
struct kdleaf *createLeaf(int axis, double **points){
    struct kdleaf *kdstruct = (struct kdleaf*)malloc(sizeof(struct kdleaf));
    kdstruct->axis = axis;
    kdstruct->origin = points;
    return kdstruct;
}

//Create full recursive structure of kdtree
struct kdleaf *CreateStructure(double **points,struct kdleaf **headkdstruct,int axis,int *length,int lower,int upper){//,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    if(axis == *length){
        axis = 0;//will correct back to right axis
    }
    struct kdleaf *currentkdstruct = *headkdstruct;
    if(upper-lower > 5){
        //Too many points, we will subdivide
        //first sorts points
        skeleQuickSort(points,lower,upper,axis,length);
        //choses node
        int nodeindex = (upper + lower) / 2; 
        double **nodepoint = (double**)malloc(sizeof(double*));
        nodepoint[0] = points[nodeindex];
        //make leaf & allocate
        currentkdstruct = createLeaf(axis,nodepoint);
        bool flag = false;
        currentkdstruct->flag = flag; 
        
        //points now sorted, Next applying our node point and split the list to the appropiate list

        currentkdstruct->left = CreateStructure(points,&currentkdstruct->left,axis + 1,length,lower,nodeindex - 1);
        currentkdstruct->right = CreateStructure(points,&currentkdstruct->right,axis+ 1,length,nodeindex + 1,upper);
        *headkdstruct = currentkdstruct;
    }
    else if(upper-lower > 0){
        //big enough for a layer, creates hold for collection of points
        double **hold = (double**)malloc((upper - lower) * sizeof(double*));
        for(int i = 0; i < upper - lower; i++){
            hold[i] = points[lower + i];
        }
        //makes leaf
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
    //we save each point when it gets finished to help prevent memory error
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

void makeSkeleton(double **points,struct kdleaf *kdstruct,int *dim,int *length,double *mindis,char path[80],bool isvofactive){
    int MAXCYCLES = 50;

    //allocate needed space
    double guessr = *length;
    double **skeleton = (double**)malloc((*length + 1) * sizeof(double*));
    for(int i = 0; i < *length+1;i++){
        skeleton[i] = (double*)malloc((*dim + 1)*sizeof(double));
        for(int j = 0; j < *dim + 1; j ++){
            skeleton[i][j] = 0.;
        }
    }
    double **centerPoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    double *radius = (double*)malloc(MAXCYCLES * sizeof(double));
    double **interfacePoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    for(int i = 0; i < *length + 1; i++){
        //Goes through each point of the list & generates a skeleton point
        //First step is to make an initial guess
        bool completeCase = true;
        int index = 0;
        
        //we Grab our temp centerpoint and our ignore point
        double *ttpoint = (double*)malloc(*dim * sizeof(double));
        double **ilist = (double**)malloc(sizeof(double*));
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
        }
        else{
            double x = points[i][0] - points[i][3] * guessr;
            double y = points[i][1] - points[i][4] * guessr;
            double z = points[i][2] - points[i][5] * guessr;
            double ix = points[i][0];
            double iy = points[i][1];
            double iz = points[i][2];
            ttpoint[0] = x;
            ttpoint[1] = y;
            ttpoint[2] = z;
            ignorepoint[0] = ix;
            ignorepoint[1] = iy;
            ignorepoint[2] = iz;
        }
        ilist[0] = ignorepoint;
        int ileng = 1;
        double lowestdistance = 0; 
        centerPoint[index] = ttpoint;//this gets overwritten first step, not sure if this plays a role, but this centerpoint will always be completely wrong
        
        //We calculate our starting furthest interface point
        interfacePoint[index] = getNearest(ttpoint,kdstruct,length,dim,ilist,&ileng,&lowestdistance);
        
        //find starting radius for our starting interface point
        double *sendpoint = points[i];
        radius[index] = getRadius(sendpoint,interfacePoint[index],dim);
        
        //now we itterate until convergence
        while(completeCase){
            //get timestep centerpoint and ignore point(not completely nessicary)
            if(*dim == 2){
                double x = points[i][0] - points[i][2] * radius[index];
                double y = points[i][1] - points[i][3] * radius[index];
                double ix = points[i][0];
                double iy = points[i][1];
                ttpoint[0] = x;
                ttpoint[1] = y;
                ignorepoint[0] = ix;
                ignorepoint[1] = iy;
            }
            else{
                double x = points[i][0] - points[i][3] * radius[index];
                double y = points[i][1] - points[i][4] * radius[index];
                double z = points[i][2] - points[i][5] * radius[index];
                double ix = points[i][0];
                double iy = points[i][1];
                double iz = points[i][2];
                ttpoint[0] = x;
                ttpoint[1] = y;
                ttpoint[2] = z;
                ignorepoint[0] = ix;
                ignorepoint[1] = iy;
                ignorepoint[2] = iz;
            }
            ilist[0] = ignorepoint;
            int ileng = 1;
            //calculate our centerpoint
            centerPoint[index] = ttpoint;
            
            //calculate our interface point closest to the last centerpoint
            lowestdistance = 0;
            interfacePoint[index + 1] = getNearest(centerPoint[index],kdstruct,length,dim,ilist,&ileng,&lowestdistance);
            
            //finds the radius of our point and interface point 
            radius[index + 1] = getRadius(points[i],interfacePoint[index + 1],dim);
            
            //get distance comp, abs distance from point->interface point, for converge check
            double distancecomp = getDistance(interfacePoint[index + 1],points[i],dim);
            
            //check for completion of skeleton point
            if(radius[index] != 0. && fabs(radius[index] - radius[index + 1]) < *mindis){
                //convergance conditions
                //our center point should remain the same
                for(int ii = 0; ii < *dim;ii++){
                    skeleton[i][ii] = centerPoint[index][ii];
                }
                skeleton[i][*dim] = radius[index + 1];
                outputskeleton(skeleton[i],dim,path);
                completeCase = false;
            }
            //temp out for smooth 
            if(isvofactive){
                if(index > 0 && distancecomp < radius[index + 1]){
                    //distance of point->interface point is less than our radius, so we want to backstep
                    for(int ii = 0; ii < *dim;ii++){
                        skeleton[i][ii] = centerPoint[index-1][ii];
                    }
                    skeleton[i][*dim] = radius[index];
                    outputskeleton(skeleton[i],dim,path);
                    completeCase = false;
                }
                else if(radius[index + 1] < *mindis){
                    //distance of point->interface point is less than our radius, so we want to backstep
                    for(int ii = 0; ii < *dim;ii++){
                        skeleton[i][ii] = centerPoint[index-1][ii];
                    }
                    //skeleton[i] = centerPoint[index-1];
                    skeleton[i][*dim] = radius[index];
                    outputskeleton(skeleton[i],dim,path);
                    completeCase = false;
                }

            }
            index += 1;
            if(index >= MAXCYCLES){
                fprintf(stdout,"broken\n");
                break;
            }
        }
        free(ttpoint);
        free(ignorepoint);
        free(ilist);
    }
    //free up needed values to prevent memory error
    free(radius);
    free(centerPoint);
    free(interfacePoint);
    for(int i = 0; i < *length+1;i++){
        free(skeleton[i]);
    }
    free(skeleton);
    radius = NULL;
    centerPoint = NULL;
    interfacePoint = NULL;
    skeleton = NULL;
}
void skeletize(double **points,int *length,int *dim,char path[80],double *mindis,bool isvofactive){
    
    //Create our kd tree for calculation
    struct kdleaf *kdstruct;
    CreateStructure(points,&kdstruct,0,dim,0,*length);//make kd-struct
    
    //Next our skeleton will be calculated
    makeSkeleton(points,kdstruct,dim,length,mindis,path,isvofactive);
    
    //Finally clean up to prevent memeory error
    for(int i = 0;i < *length + 1; i++){
        free(points[i]);
    }
    free(points);
    points = NULL;
    kdDestroy(kdstruct);
}
