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
    //fprintf(stdout,"p1x=%f p1y=%f\n",point1[0],point1[1]);
    //fprintf(stdout,"p2x=%f p2y=%f\n",point2[0],point2[1]);
    for(int i = 0; i < *dim; i++){
        double tsep = point1[i] - point2[i];
        double tabs = fabs(tsep);
        double tpow = pow(tabs,2);
        distance += tpow;
        //fprintf(stdout,"inside=%d\n",abs((int)(point1[i] - point2[i])));
        //fprintf(stdout,"p1 %f - p2 %f\n",point1[i],point2[i]);
        //fprintf(stdout,"distance=%f\n",distance);
    }
    distance = sqrt(distance);
    //fprintf(stdout,"sqrt(distance)=%f\n",distance);
    return distance;
}

double *getNearest(double *searchpoint,struct kdleaf *kdtree, int *length,int *dim,double *ignorepoint){
    double *retpoint;
    if(*dim == 2){
        retpoint = (double*)malloc(3 * sizeof(double));
    }
    else{
        retpoint = (double*)malloc(4 * sizeof(double));
    }
    if(kdtree->flag){
        //Lowest search level.
        double lowdis = 0.0;
        //fprintf(stdout,"length: %d\n",kdtree->leng);
        for(int i = 0; i < kdtree->leng; i++){
            //fprintf(stdout,"BruteForce:\n");
            //fprintf(stdout,"x:%f y=%f\n",kdtree->origin[i][0],kdtree->origin[i][1]);
            //fprintf(stdout,"ignore(%f,%f)\n",ignorepoint[0],ignorepoint[1]);
            //fprintf(stdout,"test(%f,%f)\n",kdtree->origin[i][0],kdtree->origin[i][1]);
            double tdis = getDistance(searchpoint,kdtree->origin[i],dim);
            if((lowdis == 0.0 || tdis < lowdis) && getDistance(ignorepoint,kdtree->origin[i],dim) != 0.0){
                lowdis = tdis;
                retpoint[0] = kdtree->origin[i][0];
                retpoint[1] = kdtree->origin[i][1];
                if(*dim == 3){
                    retpoint[2] = kdtree->origin[i][2];
                    retpoint[3] = lowdis;
                } 
                else{
                    retpoint[2] = lowdis;
                }
            }
        }
        //fprintf(stdout,"final dis = %f\n",lowdis);
    }
    else{
        //Can go deeper, finds node to dive to
        //fprintf(stdout,"axis = %d\n",kdtree->axis);
        //fprintf(stdout,"\n");
        //fprintf(stdout,"searchpoint x:%f y:%f\n",searchpoint[0],searchpoint[1]);
        //fprintf(stdout,"treepoint x:%f y:%f\n",*kdtree->origin[0],*kdtree->origin[1]);
        //fprintf(stdout,"axis: %d\n",kdtree->axis);
        //fprintf(stdout,"\n");
        bool side = searchpoint[kdtree->axis] < kdtree->origin[0][kdtree->axis] ;
        if(side){
            retpoint = getNearest(searchpoint,kdtree->left,length,dim,ignorepoint);
            //fprintf(stdout,"first left got: %f\n",retpoint[2]);
        }
        else{
            retpoint = getNearest(searchpoint,kdtree->right,length,dim,ignorepoint);   
            //fprintf(stdout,"first right got: %f\n",retpoint[2]);
        }
        //next we will test if we need to go to the other side of the tree
        double axisnodedis = abs((int)(searchpoint[kdtree->axis] - kdtree->origin[0][kdtree->axis]));//This only measure along one axis
        if(axisnodedis < retpoint[*dim]){
            double *tempretpoint;
            if(side){
                tempretpoint = getNearest(searchpoint,kdtree->right,length,dim,ignorepoint);
                //fprintf(stdout,"second right got: %f\n",tempretpoint[2]);
            }
            else{
                tempretpoint = getNearest(searchpoint,kdtree->left,length,dim,ignorepoint);   
                //fprintf(stdout,"second left got: %f\n",tempretpoint[2]);
            }
            if(tempretpoint[*dim] < retpoint[*dim]){
                //fprintf(stdout,"swap\n");
                retpoint = tempretpoint;
            }
        }
        //finally check the node
        //fprintf(stdout,"NodeCheck:\n");
        //fprintf(stdout,"x:%f y=%f\n",*kdtree->origin[0],*kdtree->origin[1]);
        double nodedis = getDistance(searchpoint,kdtree->origin[0],dim);
        //fprintf(stdout,"node dis = %f\n",nodedis);
        //fprintf(stdout,"\n");
        if(getDistance(ignorepoint,kdtree->origin[0],dim) != 0.0 && nodedis < retpoint[*dim]){
            retpoint[0] = kdtree->origin[0][0];
            retpoint[1] = kdtree->origin[0][1];
            if(*dim == 3){
                retpoint[2] = kdtree->origin[0][2];
                retpoint[3] = nodedis;
            } 
            else{
                retpoint[2] = nodedis;
            }
        }
    }
    //fprintf(stdout,"final got: %f\n",retpoint[2]);
    return retpoint;
}


//Destroy
void kdDestroy(struct kdleaf *kdtree){
    //not working
    fprintf(stdout,"destroying...\n");
    struct kdleaf kdtree2 = *kdtree;
    free(kdtree2.origin);
    if(kdtree2.flag){
        fprintf(stdout,"destroyed bot\n");
    }
    else{
        fprintf(stdout,"destroyed mid\n");
        kdDestroy(kdtree2.left);
        free(kdtree2.left);
        kdDestroy(kdtree2.right);
        free(kdtree2.right);
    }
}
struct kdleaf *createLeaf(int axis, double **points){
    struct kdleaf *kdtree = (struct kdleaf*)malloc(sizeof(struct kdleaf));
    kdtree->axis = axis;
    kdtree->origin = points;
    return kdtree;
}
struct kdleaf *CreateStructure(double **points,struct kdleaf **headkdtree,int axis,int *length,int lower,int upper){//,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    //tree.leftbound = leftbound;
    //tree.rightbound = rightbound;
    if(axis == *length){
        axis = 0;//will correct back to right axis
    }
    struct kdleaf *currentkdtree = *headkdtree;
    if(upper-lower > 5){
        //fprintf(stdout,"size of numbers = %d\n",*length);
        
        skeleQuickSort(points,lower,upper,axis,length);
        int nodeindex = (upper + lower) / 2; 
        //fprintf(stdout,"main create %d %d\n",upper,lower);
        double **nodepoint = (double**)malloc(sizeof(double*));
        nodepoint[0] = points[nodeindex];
        currentkdtree = createLeaf(axis,nodepoint);
        currentkdtree->flag = false; 
        
        //points now sorted, Next applying our node point and split the list to the appropiate list
        //getting new bounds for lower layers
        //bounds change depending on axis dividing

        //fprintf(stdout,"making new left: %d %d \n",lower,nodeindex-1);
        currentkdtree->left = CreateStructure(points,&currentkdtree->left,axis + 1,length,lower,nodeindex - 1);
        //fprintf(stdout,"making new right: %d %d \n",nodeindex+1,upper);
        currentkdtree->right = CreateStructure(points,&currentkdtree->right,axis+ 1,length,nodeindex + 1,upper);
        //fprintf(stdout,"main create good\n");
        *headkdtree = currentkdtree;
    }
    else if(upper-lower > 0){
        //fprintf(stdout,"alt create\n");
        double **hold = (double**)malloc((upper - lower) * sizeof(double*));
        for(int i = 0; i < upper - lower; i++){
            hold[i] = (double*)malloc(sizeof(double)*2*(*length));
            //Now we have allocated space for this
            memcpy(hold[i],points[lower + i],sizeof(double) * 2 * (*length));
        }
        currentkdtree = createLeaf(axis,hold);
        currentkdtree->flag = true;
        int leng = upper-lower;
        currentkdtree->leng = leng;
        //fprintf(stdout,"u/l = %d %d\n",upper,lower);

        //double **sliced = slice(*points,lower,upper,*length);
        //fprintf(stdout,"alt create good\n");
        
    }
    return currentkdtree;
}
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
double getRadius(double *point,double *interface,int *dim){
    //point has norm data attached
    //dim is size of data(makes life easier)
    fprintf(stdout,"\n starting radius!!\n");
    double top1 = 0.0;
    double bot = getDistance(interface,point,dim);
    //fprintf(stdout,"NORMY2:%f\n",point[3]);
    fprintf(stdout,"Point x:%f y:%f\n",point[0],point[1]);
    fprintf(stdout,"Inter x:%f y:%f\n",interface[0],interface[1]);
    fprintf(stdout,"distance=%f\n",bot);
    for(int i = 0; i < *dim; i++){
        //goes though each dimension
        double diff1 = point[i] - interface[i];
        double add1 = point[*dim+i] * diff1;
        top1 += add1;
        fprintf(stdout,"top=%f\n",top1);
    }
    double inside1 = top1/bot; 
    //fprintf(stdout,"\n");
    if(inside1 > 1){
        fprintf(stdout,"in-> 1\n");
        inside1 = 1;
    }
    if(inside1 < -1){
        fprintf(stdout,"in-> -1\n");
        inside1 = -1;
    }
    double theta1 = acos(inside1);
    double ret1 = (bot / (2 * cos(theta1)));
    double theta2 = asin(inside1);
    double ret2 = (bot / (2 * cos(theta2)));
    fprintf(stdout,"thetac=%f\n",theta1);
    fprintf(stdout,"thetas=%f\n",theta2);
    fprintf(stdout,"rad1=%f\n",ret1);
    fprintf(stdout,"rad2=%f\n",ret2);
    //fprintf(stdout,"\n");
    //fprintf(stdout,"\n");
    if(ret1 > 0 && ret1 > ret2){
        return ret1;
    }
    else if(ret2 > 0 && ret2 > ret1){
        return ret2;
    }
    else{
        fprintf(stdout,"ERROR NO RAD!\n");
        return 0.;
    }
}
double **makeSkeleton(double **points,struct kdleaf *kdtree,int *dim,int *length){
    int MAXCYCLES = 50;

    fprintf(stdout,"starting process\n");
    double guessr = *length;
    double **skeleton = (double**)malloc(*length * sizeof(double*));
    double **centerPoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    double *radius = (double*)malloc(MAXCYCLES * sizeof(double));
    double **interfacePoint = (double**)malloc(MAXCYCLES * sizeof(double*));
    for(int t = 0; t < *length; t++){
        skeleton[t] = (double*)malloc((*dim + 1) * sizeof(double));//set size of points + r
    }
    for(int t = 0; t < MAXCYCLES; t++){
        centerPoint[t] = (double*)malloc((*dim) * sizeof(double));//set size of points
        interfacePoint[t] = (double*)malloc((*dim) * sizeof(double));//set size of points
    }
    for(int i = 0; i < *length; i++){
        fprintf(stdout,"%d\n",i);
        //Goes through each element of the list & generates a skeleton point
        //First step is to make an initial guess
        bool completeCase = true;
        int index = 0;
        double *tpoint;
        double x;
        double y;
        double z;
        if(*dim == 2){
            double x = points[i][0] + guessr;
            double y = points[i][1] + guessr;
            double *ttpoint = (double*)malloc(2 * sizeof(double));
            ttpoint[0] = x;
            ttpoint[1] = y;
            tpoint = ttpoint;
            //fprintf(stdout,"x = %f or %f\n",x,tpoint[0]);
            //fprintf(stdout,"y = %f or %f\n",y,tpoint[1]);
        }
        else{
            double x = points[i][0] + guessr;
            double y = points[i][1] + guessr;
            double z = points[i][2] + guessr;
            double *ttpoint = (double*)malloc(3 * sizeof(double));
            ttpoint[0] = x;
            ttpoint[1] = y;
            ttpoint[2] = z;
            tpoint = ttpoint;
        }
        //fprintf(stdout,"x = %f \n",tpoint[0]);
        //fprintf(stdout,"y = %f \n",tpoint[1]);
        //We get a interface point and calulate the raidus to begin
        double *ignorepoint = (double*)malloc(*dim * sizeof(double));
        ignorepoint[0] = points[i][0];
        ignorepoint[1] = points[i][1];
        if(*dim == 3){
            ignorepoint[2] = points[i][2];
        }
        interfacePoint[index] = getNearest(tpoint,kdtree,length,dim,ignorepoint);
        //fprintf(stdout,"center= x:%f y:%f \n",tpoint[0],tpoint[1]);
        //fprintf(stdout,"skeleton for point%d x:%f y:%f nx:%f ny:%f\n",i,points[i][0],points[i][1],points[i][2],points[i][3]);
        //fprintf(stdout,"intpoint%d x:%f y:%f\n",index,interfacePoint[index][0],interfacePoint[index][1]); 
        double* sendpoint = points[i];
        //fprintf(stdout,"NORMY:%f\n",sendpoint[3]);
        radius[index] = getRadius(sendpoint,interfacePoint[index],dim);
        while(completeCase){
            //first calculate new center point based on the last radius given
            double *tpoint;
            if(*dim == 2){
                x = points[i][0] - radius[index] * points[i][2];
                y = points[i][1] - radius[index] * points[i][3];
                double *ttpoint = (double*)malloc(2 * sizeof(double));
                ttpoint[0] = x;
                ttpoint[1] = y;
                tpoint = ttpoint;
                //fprintf(stdout,"x = %f \n",x);
                //fprintf(stdout,"y = %f \n",y);
            }
            else{
                x = points[i][0] - radius[index] * points[i][3];
                y = points[i][1] - radius[index] * points[i][4];
                z = points[i][2] - radius[index] * points[i][5];
                double *ttpoint = (double*)malloc(3 * sizeof(double));
                ttpoint[0] = x;
                ttpoint[1] = y;
                ttpoint[2] = z;
                tpoint = ttpoint;
            }
            centerPoint[index] = tpoint;
            //then calculate the closest point
            //fprintf(stdout,"\n");
            //fprintf(stdout,"center= x:%f y:%f \n",tpoint[0],tpoint[1]);
            //fprintf(stdout,"y1 = %f \n",tpoint[1]);
            //fprintf(stdout,"\n");
            double *ignorepoint = (double*)malloc(*dim * sizeof(double));
            ignorepoint[0] = points[i][0];
            ignorepoint[1] = points[i][1];
            if(*dim == 3){
                ignorepoint[2] = points[i][2];
            }
            interfacePoint[index + 1] = getNearest(tpoint,kdtree,length,dim,ignorepoint);
            radius[index + 1] = getRadius(points[i],interfacePoint[index + 1],dim);
            //fprintf(stdout,"rad is %f\n",radius[index + 1]);
            //fprintf(stdout,"\n");
            //fprintf(stdout,"rad (%d,%d) = %f\n",i,index,radius[index+1]);
            if(radius[index] != 0. && fabs(radius[index] - radius[index + 1]) < 0.001){
                //convergance conditions
                //our center point should remain the same
                for(int j = 0; j < *dim; j++){
                    //Create skeleton point
                    //fprintf(stdout,"center:%f\n",centerPoint[index][j]);
                    skeleton[i][j] = centerPoint[index][j];
                }
                skeleton[i][*length] = radius[index];
                fprintf(stdout,"final rad = %f\n",skeleton[i][*length] = radius[index]);
                completeCase = false;
            }
            index += 1;
            if(index >= MAXCYCLES){
                fprintf(stdout,"broken\n");
                break;
            }
        }
        //fprintf(stdout,"\n");
    }
    return skeleton;
}
void outputskeleton(double **points, int *length, int *dim, char path[80]){
    FILE *fp = fopen(path,"w");
    for(int i = 0; i < *length; i++){
        if(*dim == 2){
            //fprintf(stdout,"x:%f\n",points[i][0]);
            //fprintf(stdout,"y:%f\n",points[i][1]);
            //fprintf(stdout,"r:%f\n",points[i][2]);
            fprintf(fp,"%f %f %f\n",points[i][0],points[i][1],points[i][2]);
        }
        else{
            fprintf(fp,"%f %f %f %f\n",points[i][0],points[i][1],points[i][2],points[i][3]);
        }
    }
    fclose(fp);
}
void skeletize(double **points,int *length,int *dim,char path[80]){
    //initial setup
    struct kdleaf *kdtree;
    //sizes
    //int length = (*(&points[0] + 1) - points[0]) / 2;//gets dimension of skeleton data, should give 2 or 3
    //fprintf(stdout,"length(dim) is: %d\n",*length);
    //int max = *(&points + 1) - points;//gets amount of points being put in
    //fprintf(stdout,"we have %d points\n",*max);
    CreateStructure(points,&kdtree,0,dim,0,*length);//make kd-tree
    fprintf(stdout,"we have completed kdtree\n");
    //Next we will skeletonize all the points in our list
    double **skeleton = makeSkeleton(points,kdtree,dim,length);
    outputskeleton(skeleton,length,dim,path);
    //print results
    //for(int i = 0; i < *length; i++){
    //    if(*dim == 2){
    //        fprintf(stdout,"Skeleton %d: [%f,%f] Radius: %f\n",i,skeleton[i][0],skeleton[i][1],skeleton[i][2]); 
    //    }
    //    else{
    //        fprintf(stdout,"Skeleton %d: [%f,%f,%f] Radius: %f\n",i,skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]); 
    //    }
    //}
    //kdDestroy(kdtree);
}

