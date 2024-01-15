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
    bool *flag;
    int *axis;
    int *leng;
    double **origin;//node point [x,y] or [x,y,z] @ pos 0
                    //If lowest layer, then list of points
    struct kdleaf *left,*right; 
};

//method for getting distance
double getDistance(double *point1,double *point2){
    //sqrt sum of squares for distance
    double distance = 0;
    for(int i = 0; i < dimension; i++){
        double tsep = fabs(point1[i] - point2[i]);
        double tpow = pow(tsep,2);
        distance += tpow;
    }
    distance = sqrt(distance);
    return distance;
}

double *getNearest(double *searchpoint,struct kdleaf *kdstruct, int *length,double **ignorepoint,int *ileng,double *lowestdistance){
    double *retpoint = NULL;
    if(kdstruct->flag == NULL){
        printf("error\n");
    }
    else{
        if(kdstruct->leng == NULL){
            printf("error error\n");
        }
    }
    if(*kdstruct->flag && *kdstruct->leng != 0){
        //Lowest search level.
        double lowdis = 0.0;
        for(int i = 0; i < *kdstruct->leng; i++){
            double tdis = getDistance(searchpoint,(kdstruct->origin)[i]);
            if(*ileng == 1){
                if((lowdis == 0.0 || tdis < lowdis) && getDistance(ignorepoint[0],(kdstruct->origin)[i]) != 0.0){
                    lowdis = tdis;
                    retpoint = (kdstruct->origin)[i];
                    *lowestdistance = lowdis;
                }
            }
            else{
                bool passes = true;
                for(int j = 0; j < *ileng; j++){
                    if((*kdstruct->origin)[i] == *ignorepoint[j]){
                        passes = false;
                        break;
                    }
                }
                if((lowdis == 0.0 || tdis < lowdis) && passes){
                    lowdis = tdis;
                    retpoint = (kdstruct->origin)[i];
                    *lowestdistance = lowdis;
                }
            }
        }
    }
    else if(!(*kdstruct->flag)){
        //can go deeper
        bool side = searchpoint[*kdstruct->axis] < (kdstruct->origin)[0][*kdstruct->axis];
        if(side){
            if(kdstruct->left != NULL){
                retpoint = getNearest(searchpoint,kdstruct->left,length,ignorepoint,ileng,lowestdistance);
            }
        }
        else{
            if(kdstruct->right != NULL){
                retpoint = getNearest(searchpoint,kdstruct->right,length,ignorepoint,ileng,lowestdistance);   
            }
        }
        //next we will test if we need to go to the other side of the struct
        double axisnodedis = fabs(searchpoint[*kdstruct->axis] - (kdstruct->origin)[0][*kdstruct->axis]);//This only measure along one axis
        if(axisnodedis < *lowestdistance || *lowestdistance == 0.0){
            double *tempretpoint;
            double templowdis = 0.0;
            if(side){
                if(kdstruct->right != NULL){
                    tempretpoint = getNearest(searchpoint,kdstruct->right,length,ignorepoint,ileng,&templowdis);   
                }
            }
            else{
                if(kdstruct->left != NULL){
                    tempretpoint = getNearest(searchpoint,kdstruct->left,length,ignorepoint,ileng,&templowdis);
                }
            }
            if((*lowestdistance > templowdis|| *lowestdistance == 0.0) && templowdis != 0.0){
                retpoint = tempretpoint;
                *lowestdistance = templowdis;
            }
        }
        //finally check the node
        double nodedis = getDistance(searchpoint,(kdstruct->origin)[0]);
        if(*ileng == 1){
            if(getDistance(ignorepoint[0],(kdstruct->origin)[0]) != 0.0 && (nodedis < *lowestdistance || *lowestdistance == 0.0)){
                retpoint = (kdstruct->origin)[0];
                *lowestdistance = nodedis;
            }
        }
        else{
            bool passes = true;
            for(int j = 0; j < *ileng; j++){
                if((*kdstruct->origin)[0] == *ignorepoint[j]){
                    passes = false;
                    break;
                }
            }
            if((nodedis < *lowestdistance || *lowestdistance == 0.0) && passes){
                retpoint = (kdstruct->origin)[0];
                *lowestdistance = nodedis;
            }
        }
    }
    return retpoint;
}

//Destroy tree
void kdDestroy(struct kdleaf **kdstruct){
    if(*kdstruct != NULL){
        struct kdleaf *tstruct = *kdstruct;
        if(!(*(tstruct->flag))){
            if((*kdstruct)->left != NULL){
                kdDestroy(&(*kdstruct)->left);
            }
            if((*kdstruct)->right != NULL){
                kdDestroy(&(*kdstruct)->right);
            }
        }
        if((*kdstruct)->origin != NULL){
            for(int i = 0; i < *(tstruct->leng); i++){
                free((*kdstruct)->origin[i]);
            }
            free((*kdstruct)->origin);
        }
        free((*kdstruct)->leng);
        free((*kdstruct)->axis);
        free((*kdstruct)->flag);
        free(*kdstruct);
        *kdstruct = NULL;
    }
}

//Create full recursive structure of kdtree
void CreateStructure(double **points,struct kdleaf **kdstruct,int axis,int lower,int upper,int editleng){//,double *leftbound, double *rightbound){
    //First we will sort along each needed axis
    struct kdleaf* headkdstruct = malloc(sizeof(struct kdleaf));
    if(axis == dimension){
        axis = 0;//will correct back to right axis
    }
    //struct kdleaf currentkdstruct = *headkdstruct;
    if(upper-lower > 5){
        //Too many points, we will subdivide
        //first sorts points
        skeleQuickSort(points,lower,upper,axis,editleng);
        //choses node
        int nodeindex = (upper + lower) / 2; 
        headkdstruct->origin = malloc(sizeof(double*));
        headkdstruct->origin[0] = malloc(sizeof(double) * (dimension + editleng));
        if(headkdstruct->origin == NULL){
            printf("errL1\n");
        }
        headkdstruct->origin[0][0] = points[nodeindex][0];
        headkdstruct->origin[0][1] = points[nodeindex][1];
        if(dimension + editleng == 3){
            headkdstruct->origin[0][2] = points[nodeindex][2];
        }
        if(headkdstruct->origin[0] == NULL){
            printf("errL2\n");
        }
        //make leaf & allocate
        int *taxis = malloc(sizeof(int));
        *taxis = axis;
        headkdstruct->axis = taxis;
        bool *flag = malloc(sizeof(bool));
        *flag = false;
        headkdstruct->flag = flag; 
        headkdstruct->leng = malloc(sizeof(int));
        *headkdstruct->leng = 1;
        //points now sorted, Next applying our node point and split the list to the appropiate list

        CreateStructure(points,&headkdstruct->left,axis + 1,lower,nodeindex - 1,editleng);
        CreateStructure(points,&headkdstruct->right,axis+ 1,nodeindex + 1,upper,editleng);
    }
    else{
        if (upper - lower > 0){
            headkdstruct->origin = malloc(sizeof(double*) * (upper - lower) );
            if(headkdstruct->origin == NULL){
                printf("errL1\n");
            }
            for(int i = 0; i < upper - lower; i++){
               
                headkdstruct->origin[i] = malloc(sizeof(double) * (dimension + editleng));     
                if(headkdstruct->origin[i] == NULL){
                    printf("errL2\n");
                }
                headkdstruct->origin[i][0] = points[lower + i][0];
                headkdstruct->origin[i][1] = points[lower + i][1];
#if dimension == 3    
                headkdstruct->origin[i][2] = points[lower + i][2];
#endif
            }
            //makes leaf
            int *taxis = malloc(sizeof(int));
            *taxis = axis;
            headkdstruct->axis = taxis;
            bool *flag = malloc(sizeof(bool));
            *flag = true;
            headkdstruct->flag = flag;
            int leng = upper-lower;
            headkdstruct->leng = malloc(sizeof(int));
            *headkdstruct->leng = leng;
        }
        else{
            printf("error c3\n");
            headkdstruct = NULL;
            return;
        }
    }
    *kdstruct = headkdstruct;
}

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

double getRadius(double *point,double *interface){
    //point has norm data attached
    //dim is size of data(makes life easier)
    double top1 = 0.0;
    double bot = getDistance(interface,point);
    for(int i = 0; i < dimension; i++){
        //goes though each dimension
        double diff1 = point[i] - interface[i];
        double add1 = point[dimension+i] * diff1;
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
        return -ret1;
    }
}

void outputskeleton(double *points, double *interf,double alpha , int indx,char path[80]){
    //we save each point when it gets finished to help prevent memory error
    FILE *fp = fopen(path,"a");
    //output with alpha
    if(dimension == 2){
        fprintf(fp,"%f %f %f %f %d\n",points[0],points[1],points[2],alpha,indx);
    }
    else{
        fprintf(fp,"%f %f %f %f %f %d\n",points[0],points[1],points[2],points[3],alpha,indx);
    }
    fflush(fp);
    fclose(fp);
}

void makeSkeleton(double **points,struct kdleaf *kdstruct,int *length,double *mindis,char path[80],bool isvofactive,double *disRatio,int *newl,double ***pskeleton){
    int MAXCYCLES = 50;
    //allocate needed space
    double guessr = *length;
    int extra = 3;//save r, alpha, kappa to relevant skeletonpoints
    double **skeleton = malloc((*length + 1) * sizeof(double*));
    for(int i = 0; i < *length+1;i++){
        skeleton[i] = calloc((dimension + extra) , sizeof(double));
    }
    double **centerPoint = malloc(MAXCYCLES * sizeof(double*));
    double *radius = malloc(MAXCYCLES * sizeof(double));
    double **interfacePoint = malloc(MAXCYCLES * sizeof(double*));
    int captured = 0;
    for(int i = 0; i < *length + 1; i++){
        //printf("skel %d / %d\n",i,*length);
        //Goes through each point of the list & generates a skeleton point
        //First step is to make an initial guess
        bool completeCase = true;
        int index = 0; 
        //we Grab our temp centerpoint and our ignore point
        double *ttpoint = malloc(dimension * sizeof(double));
        double **ilist = malloc(sizeof(double*));
        double *ignorepoint = malloc(dimension * sizeof(double));
        double ix = points[i][0];
        double iy = points[i][1];
        ignorepoint[0] = ix;
        ignorepoint[1] = iy;
#if dimension == 2    
        double x = points[i][0] - points[i][2] * guessr;
        double y = points[i][1] - points[i][3] * guessr;
        ttpoint[0] = x;
        ttpoint[1] = y;
#else    
        double x = points[i][0] - points[i][3] * guessr;
        double y = points[i][1] - points[i][4] * guessr;
        double z = points[i][2] - points[i][5] * guessr;
        double iz = points[i][2];
        ttpoint[0] = x;
        ttpoint[1] = y;
        ttpoint[2] = z;
        ignorepoint[2] = iz;
#endif
        ilist[0] = ignorepoint;
        int ileng = 1;
        double lowestdistance = 0; 
        centerPoint[index] = ttpoint;//this gets overwritten first step, not sure if this plays a role, but this centerpoint will always be completely wrong 
        //We calculate our starting furthest interface point
        interfacePoint[index] = getNearest(ttpoint,kdstruct,length,ilist,&ileng,&lowestdistance);
        //find starting radius for our starting interface point
        double *sendpoint = points[i];
        radius[index] = getRadius(sendpoint,interfacePoint[index]); 
        //now we itterate until convergence
        while(completeCase){
            if(index > MAXCYCLES - 2){
                break;
            }
            //get timestep centerpoint and ignore point(not completely nessicary)
            if(dimension == 2){
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
            interfacePoint[index + 1] = getNearest(centerPoint[index],kdstruct,length,ilist,&ileng,&lowestdistance);
            //finds the radius of our point and interface point 
            radius[index + 1] = getRadius(points[i],interfacePoint[index + 1]);
            //get distance comp, abs distance from point->interface point, for converge check
            double distancecomp = getDistance(interfacePoint[index + 1],points[i]);
            double alpha = distancecomp / radius[index + 1];
            //check for completion of skeleton point
            if(radius[index] != 0. && fabs(radius[index] - radius[index + 1]) < *mindis){
                //convergance conditions
                //our center point should remain the same
                for(int ii = 0; ii < dimension;ii++){
                    skeleton[i][ii] = centerPoint[index][ii];
                }
                skeleton[i][dimension] = radius[index + 1];
                skeleton[i][dimension+1] = alpha;
                skeleton[i][dimension+2] = points[i][(dimension * 2)]; 
                captured++;
                completeCase = false;
            }
            //temp out for smooth 
            if(isvofactive){
                if(index > 0 && distancecomp < radius[index + 1]){
                    //distance of point->interface point is less than our radius, so we want to backstep
                    for(int ii = 0; ii < dimension;ii++){
                        skeleton[i][ii] = centerPoint[index-1][ii];
                    }
                    skeleton[i][dimension] = radius[index];
                    skeleton[i][dimension+1] = alpha;
                    skeleton[i][dimension+2] = points[i][(dimension * 2)]; 
                    captured++;
                    completeCase = false;
                }
                else if(radius[index + 1] < *mindis){
                    //distance of point->interface point is less than our radius, so we want to backstep
                    for(int ii = 0; ii < dimension;ii++){
                        skeleton[i][ii] = centerPoint[index-1][ii];
                    }
                    skeleton[i][dimension] = radius[index];
                    skeleton[i][dimension+1] = alpha;
                    skeleton[i][dimension+2] = points[i][(dimension * 2)]; 
                    captured++;
                    completeCase = false;
                }

            }
            else{
                if(radius[index + 1] < *mindis){
                    completeCase = false;
                }
            }
            index = index + 1;
        }
        free(ttpoint);
        free(ignorepoint);
        free(ilist);
    }
    //free up needed values to prevent memory error
    free(radius);
    free(centerPoint);
    free(interfacePoint);
    radius = NULL;
    centerPoint = NULL;
    interfacePoint = NULL;
    *newl = captured;
    *pskeleton = skeleton;
}
void skeletize(double **points,int *length,char path[80],double *mindis,bool isvofactive,double wantedAngle,double ***pskeleton){
    double **skeleton = NULL; 
    if(*length > 0){
        //Create our kd tree for calculation
        struct kdleaf *kdstruct = NULL;
        CreateStructure(points,&kdstruct,0,0,*length,0);//make kd-struct
        if(kdstruct== NULL){
            printf("error\n");
        }
        //Calculate our distance ratio
        if(wantedAngle > 180){
            wantedAngle = 180;
        }
        double calcratio = 0.0;
        calcratio = tan(wantedAngle);
        //Next our skeleton will be calculated
        int newl = 0;
        makeSkeleton(points,kdstruct,length,mindis,path,isvofactive,&calcratio,&newl,&skeleton);
        //Finally clean up to prevent memeory error
        for(int i = 0;i < *length + 1; i++){
            free(points[i]);
        }
        free(points);
        points = NULL;
        kdDestroy(&kdstruct);
        length = &newl;

    }
    *pskeleton = skeleton;
}
void thinSkeleton(double ***pskeleton,int *length,double *alpha,double *thindis,char outname[80]){
    double **skeleton = *pskeleton;
    //we need to handle situations 
    printf("inl = %d / %f\n",*length,*thindis);
    FILE * fp = fopen (outname, "w");
    if(*length > 0){
        int holdl = *length;
        bool addq = false;
        for(int i = *length - 1; i >=0; i--){
            //printf("%f\n",*thindis);
            if((skeleton[i][dimension+1] < *alpha) || skeleton[i][dimension] < *thindis || skeleton[i][dimension] > 0.5){// || (skeleton[i][dimension]  > *thindis) || (skeleton[i][dimension+2] < 1.)){
                //If bad point we will shift everything down one; removing it later
#if dimension == 2
                fprintf(fp,"%f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4]);
#else
                fprintf(fp,"%f %f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4],skeleton[i][5]);
#endif
                if(!addq){
                    addq = true;
                }
                for(int j = i+1; j < *length; j++){
                    skeleton[j-1][0] = skeleton[j][0];//x
                    skeleton[j-1][1] = skeleton[j][1];//y
                    skeleton[j-1][2] = skeleton[j][2];//z/r
                    skeleton[j-1][3] = skeleton[j][3];//r/a
                    skeleton[j-1][4] = skeleton[j][4];//a/k
#if dimension == 3                    
                    skeleton[j-1][5] = skeleton[j][5];//k
#endif
                }
                int L = *length - 1;
                *length = L;
            }
        }
        if(*length > 0){
            for(int i = *length; i < holdl; i++){
                free(skeleton[i]);
            }
            skeleton = realloc(skeleton,*length * sizeof(double*));
        }
        else{
            for(int i = 0; i < holdl + 1; i++){
                free(skeleton[i]);
            }
            free(skeleton);
            skeleton = NULL;
        }
    }
    fflush(fp);
    fclose(fp);
    printf("outl = %d\n",*length);
    *pskeleton = skeleton;
}
//Gets spline
void getCoeffGE(int row, int col, double ***a, double **x){
    int i,j,k;
    for(i = 0; i < row - 1; i++){
        //swap
        for(k = i + 1; k < row; k++){
            if(fabs((*a)[i][i])<fabs((*a)[k][i])){
                for(j = 0; j < col; j++){
                    double temp;
                    temp = (*a)[i][j];
                    (*a)[i][j] = (*a)[k][j];
                    (*a)[k][j] = temp;
                }
            }
        }
        //Gauss
        for(k = i + 1; k <row;k++){
            double term = (*a)[k][i] / ((*a)[i][i] + SEPS);
            for(j = 0; j < col; j++){
                (*a)[k][j] = (*a)[k][j] - term * (*a)[i][j];
            }
        }
    }
    //sub
    for(i = row - 1; i >= 0; i--){
        (*x)[i]=(*a)[i][col-1];
        for(j  = i + 1; j < col - 1; j++){
            (*x)[i] = (*x)[i] - (*a)[i][j] * (*x)[j];
        }
        (*x)[i] = (*x)[i] / ((*a)[i][i] + SEPS);
    }
}
//Spline
struct skeleBounds{
    double *x;//By convention x,y is the bottom left corner of the bounds;
    double *y;
    double *z;
    double **points;
    int *leng;
    double *density;
    bool *hasNode;
    double *nodepoint;
    double **roi;//Region of Interest [xmin,ymin,zmin],[xmax,ymax,zmax]
    int *smode;//selection mode
    int *connections;
    int *closedis1;//closest distance to node point
    int *closedis2;//closest distance to node point
    int *closeid1;//closest node's id
    int *closeid2;//closest node's id
};
struct skeleDensity{
    //overall struct, this holds relevant values needed to be calculated &
    //our grid
    struct skeleBounds ***sB;//we have a 3D matrix of skele bounds refrenced by [row][col][dep]
    //NOTE: dep is only deeper than 1 if dim = 3
    int *row;
    int *col;
    int *dep;
    int *ncount;
    int *pcount;
    double *xmax;
    double *xmin;
    double *ymax;
    double *ymin;
    double *zmax;
    double *zmin;
    double *dx;
    double *dy;
    double *dz;
};

void createSD(struct skeleDensity **sd,double **inpts,int *inleng,int *dim,double delta,double xmax,double xmin,double ymax,double ymin,double zmax,double zmin){
    //Here we allocate everything needed in our structure
    struct skeleDensity *allocsd = malloc(sizeof(struct skeleDensity));
    //x alloc & set
    double *txmax = malloc(sizeof(double)); 
    *txmax = xmax;
    allocsd->xmax = txmax;
    double *txmin = malloc(sizeof(double)); 
    *txmin = xmin;
    allocsd->xmin = txmin; 
    //y alloc & set
    double *tymax = malloc(sizeof(double)); 
    *tymax = ymax;
    allocsd->ymax = tymax; 
    double *tymin = malloc(sizeof(double)); 
    *tymin = ymin;
    allocsd->ymin = tymin; 
    //z alloc & set if needed
    if(*dim == 3){
        double *tzmax = malloc(sizeof(double)); 
        *tzmax = zmax;
        allocsd->zmax = tzmax; 
        double *tzmin = malloc(sizeof(double)); 
        *tzmin = zmin;
        allocsd->zmin = tzmin; 
    }
    int *tnc = malloc(sizeof(int));
    *tnc = 0;
    allocsd->ncount = tnc;
    int *tpc = malloc(sizeof(int));
    *tpc = 0;
    allocsd->pcount = tpc;
    //finally we will allocate & calulate our matrix of grids
    //first is dim/offset calc 
    int row = 0; int col = 0;int dep = 0;
    //fprintf(stdout,"%f-%f/%f\n",ymax,ymin,delta);
    row = (int)ceil((ymax - ymin)/delta);
    double offr = ((row * delta) - (ymax - ymin));
    col = (int)ceil((xmax - xmin)/delta);
    double offc = ((col * delta) - (xmax - xmin)); 
    double offd = 0.;
    if(*dim == 3){
        dep = (int)ceil((zmax - zmin)/delta);
        offd = ((dep * delta) - (zmax - zmin)); 
    } 
    else{
        dep = 1;
    }
    //fprintf(stdout,"offsets = i,j,k - [%f,%f,%f]\n",offr,offc,offd);
    //then we set our number of spaces
    int *trow = malloc(sizeof(int));
    int *tcol = malloc(sizeof(int));
    int *tdep = malloc(sizeof(int));
    *trow = row;
    *tcol = col;
    *tdep = dep;
    allocsd->row = trow;
    allocsd->col = tcol;
    allocsd->dep = tdep;
    
    //Next we calculate individual deltas with out new offset give the ceil
    //helps ensure equal spacings
    double *dx  = malloc(sizeof(double));
    double *dy  = malloc(sizeof(double));
    //fprintf(stdout,"x's%f,%f\n",xmax,xmin);
    *dx =  delta;
    *dy =  delta;
    allocsd->dx = dx;
    allocsd->dy = dy;
    double *dz  = malloc(sizeof(double));
    if(*dim == 3){
        *dz =  delta;
    }
    else{
        *dz =  1.;
    }
    allocsd->dz = dz;
    //fprintf(stdout,"deltas = i,j,k - [%f,%f,%f]\n",*dx,*dy,*dz);
    //then we allocate our 3D structure
    struct skeleBounds ***allocsb = malloc(row * sizeof(struct skeleBounds **));
    for(int i = 0; i < row; i++){
        allocsb[i] = malloc(col * sizeof(struct skeleBounds *));
        for(int j = 0; j < col; j++){
            allocsb[i][j] = malloc(dep * sizeof(struct skeleBounds));
        }
    }
    //Then we go through and define our x,y,z
    double xstart = xmin - offc;
    double ystart = ymin - offr;
    double zstart;
    if(*dim == 3){
        zstart = zmin  - offd;
    }
    else{
        zstart = 0.;
    }
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            for(int k = 0; k < dep; k++){
                //set up ability for having node
                bool *tb = malloc(sizeof(bool));
                *tb = false;
                allocsb[i][j][k].hasNode = tb;
                int *tsm = malloc(sizeof(int));
                *tsm = 0;
                allocsb[i][j][k].smode = tsm;
                int *tsc = malloc(sizeof(int));
                *tsc = 0;
                allocsb[i][j][k].connections = tsc;
                //allocate for distances
                int *tcd1 = malloc(sizeof(int));
                *tcd1 = -1;
                allocsb[i][j][k].closedis1 = tcd1;
                int *tcd2 = malloc(sizeof(int));
                *tcd2 = -1;
                allocsb[i][j][k].closedis2 = tcd2;
                int *tci1 = calloc(3,sizeof(int));
                allocsb[i][j][k].closeid1 = tci1;
                int *tci2 = calloc(3,sizeof(int));
                allocsb[i][j][k].closeid2 = tci2;
                //set y from row
                double *ty = malloc(sizeof (double));
                *ty = ystart + (*dy * i);
                allocsb[i][j][k].y = ty;
                //set x from col
                double *tx = malloc(sizeof (double));
                *tx = xstart + (*dx * j);
                allocsb[i][j][k].x = tx;
                //set z from dep
                double *tz = malloc(sizeof (double));
                *tz = zstart + (*dz * k);
                allocsb[i][j][k].z = tz;
                //Next Given our bounds namely .x & .x + dx ect. we sort points into the needed areas
                int *tlen = malloc(sizeof(int));
                *tlen = 0;
                double **holdpoint = malloc(*inleng * sizeof(double*));
                //fprintf(stdout,"bounds:[%f,%f,%f],[%f,%f,%f]",*tx,*ty,*tz,*tx + *dx,*ty+ *dy,*tz+ *dz);
                for(int q = 0; q < *inleng; q++){
                    //fprintf(stdout,"trying pt [%f,%f]\n",inpts[q][0],inpts[q][1]);
                    if(inpts[q][0] > *tx && inpts[q][0] <= *tx + *dx){
                        if(inpts[q][1] > *ty && inpts[q][1] <= *ty + *dy){
                            if(*dim != 3 || (inpts[q][2] > *tx && inpts[q][2] <= *tz + *dz)){
                                holdpoint[*tlen] = inpts[q];
                                *tlen = *tlen + 1;
                                //fprintf(stdout,"%d\n",*tlen);
                            }
                        }
                    }
                }
                allocsb[i][j][k].leng = tlen;//set length
                allocsb[i][j][k].points = malloc(*tlen * sizeof(double*));
                allocsb[i][j][k].roi = malloc(2 * sizeof(double*));
                (allocsb[i][j][k].roi)[0] = calloc(3 , sizeof(double));
                (allocsb[i][j][k].roi)[1] = calloc(3 , sizeof(double));
                double *calcdensity = malloc(sizeof(double));
                int extra = 1;
                if(*tlen != 0){
                    for(int q = 0; q < *dim; q++){
                        (allocsb[i][j][k].roi)[0][q] = HUGE;
                        (allocsb[i][j][k].roi)[1][q] = -HUGE;
                    }
                    for(int q = 0; q < *tlen; q++){
                        (allocsb[i][j][k].points)[q] = malloc((*dim + extra) * sizeof(double));
                        for(int p = 0; p < *dim + extra; p++){
                            (allocsb[i][j][k].points)[q][p] = holdpoint[q][p];
                            if(holdpoint[q][p] < (allocsb[i][j][k].roi)[0][p]){
                                (allocsb[i][j][k].roi)[0][p] = holdpoint[q][p];
                            }
                            if(holdpoint[q][p] > (allocsb[i][j][k].roi)[1][p]){
                                (allocsb[i][j][k].roi)[1][p] = holdpoint[q][p];
                            }
                        }
                    }
                    //adjust x
                    if((allocsb[i][j][k].roi)[0][0] - *dx / 10 > *(allocsb[i][j][k]).x){
                        (allocsb[i][j][k].roi)[0][0] = (allocsb[i][j][k].roi)[0][0] - *dx / 10;
                    }
                    else{
                        (allocsb[i][j][k].roi)[0][0] = *(allocsb[i][j][k]).x;
                    }
                    if((allocsb[i][j][k].roi)[1][0] + *dx / 10 < *(allocsb[i][j][k]).x + *dx){
                        (allocsb[i][j][k].roi)[1][0] = (allocsb[i][j][k].roi)[1][0] + *dx / 10;
                    }
                    else{
                        (allocsb[i][j][k].roi)[1][0] = *(allocsb[i][j][k]).x + *dx;
                    }
                    //adjust y
                    if((allocsb[i][j][k].roi)[0][1] - *dy / 10 > *(allocsb[i][j][k]).y){
                        (allocsb[i][j][k].roi)[0][1] = (allocsb[i][j][k].roi)[0][1] - ((*dy) / 10);
                    }
                    else{
                        (allocsb[i][j][k].roi)[0][1] = *(allocsb[i][j][k]).y;
                    }
                    if((allocsb[i][j][k].roi)[1][1] + *dy / 10 < *(allocsb[i][j][k]).y + *dy){
                        (allocsb[i][j][k].roi)[1][1] = (allocsb[i][j][k].roi)[1][1] + *dy / 10;
                    }
                    else{
                        (allocsb[i][j][k].roi)[1][1] = *(allocsb[i][j][k]).y + *dy;
                    }
                    //adjust z
                    if(*dim == 2){
                        (allocsb[i][j][k].roi)[0][2] = 0;
                        (allocsb[i][j][k].roi)[1][2] = 1;
                    }
                    else{
                        if((allocsb[i][j][k].roi)[0][2] - *dz / 10 > *(allocsb[i][j][k]).z){
                            (allocsb[i][j][k].roi)[0][2] = (allocsb[i][j][k].roi)[0][2] - *dz / 10;
                        }
                        else{
                            (allocsb[i][j][k].roi)[0][2] = *(allocsb[i][j][k]).z;
                        }
                        if((allocsb[i][j][k].roi)[1][2] + *dz / 10 < *(allocsb[i][j][k]).z + *dz){
                            (allocsb[i][j][k].roi)[1][2] = (allocsb[i][j][k].roi)[1][2] + *dz / 10;
                        }
                        else{
                            (allocsb[i][j][k].roi)[1][2] = *(allocsb[i][j][k]).z + *dz;
                        }
                    }
                    double ddx = ((allocsb[i][j][k].roi)[1][0] - (allocsb[i][j][k].roi)[0][0]);
                    double ddy = ((allocsb[i][j][k].roi)[1][1] - (allocsb[i][j][k].roi)[0][1]);
                    double ddz = ((allocsb[i][j][k].roi)[1][2] - (allocsb[i][j][k].roi)[0][2]);
                    *calcdensity = *tlen / ((ddx) * (ddy) * (ddz));//points are our 'mass' 
                    allocsb[i][j][k].density = calcdensity;
                }
                else{
                    *calcdensity = 0.;
                    allocsb[i][j][k].density = calcdensity;
                }
                //fprintf(stdout,"made?=%d\n",(allocsb[i][j][k].roi)[0] != NULL);
                free(holdpoint);//frees the extra memeory we dont actually want to keep from our sort
                //Finally we can calculate our density
                //we will correct roi to be +- delta/10 if it can
                //fprintf(stdout,"grid [%f,%f,%f] / %d\n",*dx,*dy,*dz,*tlen);
                //fprintf(stdout,"grid [%d,%d,%d] = %f\n\n",i,j,k,*calcdensity);
            }
        }
    }
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            for(int k = 0; k < dep; k++){
                if(*allocsb[i][j][k].leng > 0){
                    //fprintf(stdout,"leng1=%d\n",*(sDmain->sB[i][j][k]).leng);
                    //Here we will loop through all of our cells
                    //We need to make sure we wont scan outofbound cells
                    //So creating a check map will be beneficial
                    int len = 0;
                    for(int ic = -1; ic <= 1; ic++){
                        for(int jc = -1; jc <= 1; jc++){
                            for(int kc = -1; kc <= 1; kc++){
                                int ti = i + ic;
                                int tj = j + jc;
                                int tk = k + kc;
                                if(ti >= 0 && ti < row){
                                    if(tj >= 0 && tj < col){
                                        if(tk >= 0 && tk < dep){
                                            if(!(ic == 0 && jc == 0 && kc == 0) ){// && passcheck(ic,jc,kc,dim)){
                                                if(*(allocsb[ti][tj][tk]).leng != 0){
                                                    len++; 
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    *(allocsb[i][j][k].connections) = len;
                }
            }
        }
    }
    allocsd->sB = allocsb;
    *sd = allocsd;
    printf("complete\n");
}
void destroySD(struct skeleDensity **sd, int *dim){
    if(*sd != NULL){
        struct skeleDensity *tsd = *sd;
        //first Destroy SB's
        for(int i = 0; i < *(tsd->row);i++){
            for(int j = 0; j < *(tsd->col);j++){
                for(int k = 0; k < *(tsd->dep);k++){
                    struct skeleBounds *bounds = &(tsd->sB[i][j][k]);
                    if(*(bounds->hasNode) == true){
                        free(bounds->nodepoint);
                    }
                    free(bounds->hasNode);
                    free(bounds->x);
                    free(bounds->y);
                    free(bounds->z);
                    free(bounds->density);
                    if(bounds->points != NULL){
                        for(int q = 0; q < *(bounds->leng); q++){
                            //fprintf(stdout,"destroying %d\n",i);
                            free(bounds->points[q]);
                        }
                        free(bounds->points);
                    }
                    free(bounds->smode);
                    free(bounds->connections);
                    free(bounds->leng);
                    free(bounds->closedis1);
                    free(bounds->closedis2);
                    free(bounds->closeid1);
                    free(bounds->closeid2);
                    free(bounds->roi[0]);
                    free(bounds->roi[1]);
                    free(bounds->roi);
                }
                free(tsd->sB[i][j]);
            }
            free(tsd->sB[i]);
        }
        free(tsd->sB);
        free(tsd->row);
        free(tsd->col);
        free(tsd->dep);
        free(tsd->xmax);
        free(tsd->xmin);
        free(tsd->ymax);
        free(tsd->ymin);
        if(*dim == 3){
            free(tsd->zmax);
            free(tsd->zmin);
        }
        free(tsd->dx);
        free(tsd->dy);
        free(tsd->dz);
        free(tsd->ncount);
        free(tsd->pcount);
        free(tsd);
        *sd = NULL;
    }
    else{
        printf("error NULL\n");
    }
}

//handles deallocation of combos
void freeCombo(double ****pfindpoints,int **pcomboindex,int ***pnodeconnections,int ***pnodeindex,int *combocount,int *nicount){
    double ***findpoints = *pfindpoints;
    int **nodeindex = *pnodeindex;
    int **nodeconnections = *pnodeconnections;
    int *comboindex = *pcomboindex;
    bool clearfind = false;
    if(findpoints != NULL){
        clearfind = true;
    }
    //frees connections & points if needed
    for(int i = 0; i < *combocount; i++){
        free(nodeconnections[i]);
        if(clearfind){
            free(findpoints[i]);
        }
    }
    if(*pnodeconnections != NULL && nodeconnections != NULL){
        free(nodeconnections);
    }
    if(clearfind){
        free(findpoints);
    }
    for(int i = 0; i < *nicount; i++){
        free(nodeindex[i]);
    }
    if(nodeindex != NULL){
        free(nodeindex);
    }
    if(comboindex != NULL){
        free(comboindex);
    }
    *pfindpoints = NULL;
    *pnodeindex = NULL;
    *pnodeconnections = NULL;
    *pcomboindex = NULL;
}
void thinNodePoint(struct skeleDensity **sd,int *dim,double ****pfindpoints,int **pcomboindex,int ***pnodeconnections,int ***pnodeindex,int *combocount,int *nicount,int *mxpt, double t){
    //Thins our structure removing points with only 2 connections
    printf("\n\nThinning... watch for error\n");
    struct skeleDensity *sDmain = *sd;
    double ***findpoints = *pfindpoints;
    int **nodeindex = *pnodeindex;
    int **nodeconnections = *pnodeconnections;
    int *comboindex = *pcomboindex;
    //firstly we will allocate a count list for each index id, this will give us an idea of how  many nodes we are getting rid of and where they are located
    printf("cc=%d,ni=%d\n",*combocount,*nicount);
    for(int i = 0; i < *combocount; i++){
        printf("[%d - %d , %d ]\n",i,nodeconnections[i][0],nodeconnections[i][1]);
    }
    int *countIDs = calloc(*nicount , sizeof(int));
    int **trackIDs = malloc(*nicount * sizeof(int*));
    for(int i = 0; i < *nicount; i++){
        trackIDs[i] = calloc(4 , sizeof(int));//holds each [id00,id01,id10,id11]
    }
    int numtwo = 0;
    for(int i = 0; i < *combocount; i++){
        countIDs[nodeconnections[i][0]]++;
        countIDs[nodeconnections[i][1]]++;
    }
    //here we adjust places where we have twos touching as we done want both points, doing it before count and track allows us to stay memeory safe 
    for(int i = 0; i < *combocount; i++){
        if(countIDs[nodeconnections[i][0]] == 2 && countIDs[nodeconnections[i][1]] == 2){
            //add a fake count to position 1 in a double to prevent any redoing
            countIDs[nodeconnections[i][0]]++;
        }
    }
    for(int i = 0; i < *nicount; i++){
        if(countIDs[i] == 2){
            //Here we have a node which only has 2 connections, this means we can remove it safely without much risk of breaking the system :)
            numtwo++;
            //here we grab the useful id & assign it
            int locid = 0;
            for(int j = 0; j < *combocount; j++){
                if(nodeconnections[j][0] == i){
                    trackIDs[i][locid] = j;
                    locid++;
                    trackIDs[i][locid] = 0;
                    locid++;
                }
                else if(nodeconnections[j][1] == i){
                    trackIDs[i][locid] = j;
                    locid++;
                    trackIDs[i][locid] = 1;
                    locid++;
                }
            }
        }
    }
    printf("numtwos = %d\n",numtwo);
    //After were done counting we will alloc 
    if(numtwo > 0){
        int tnncount = 0;
        int **tempnewnode = malloc(numtwo * sizeof(int*));//for new nodeconnections
        double ***tempnewfind = malloc(numtwo * sizeof(double**));//allocs new vectors for new findpoints:)
        int *tempfindcount = calloc(numtwo, sizeof(int));//for new comboindex
        int ti = 0;
        for(int i = 0 ; i < *nicount; i++){
            if(countIDs[i] == 2){
                int countfind = (comboindex[trackIDs[i][0]]) + (comboindex[trackIDs[i][2]]);//combines count of points for proper space
                tempfindcount[ti] = countfind;
                tempnewnode[ti] = calloc(2 , sizeof(int));
                tempnewfind[ti] = malloc(countfind * sizeof(double*));
                for(int j = 0; j < countfind; j++){
                    tempnewfind[ti][j] = calloc(*dim + 1, sizeof(double));
                }
                ti++;
            }
        }
        //then we will run through each two vector and combine it into our temp
        for(int i = 0; i < *nicount; i++){
            if(countIDs[i] == 2){
                //if equal to two then we create our new values
                tempnewnode[tnncount][0] = nodeconnections[trackIDs[i][0]][!trackIDs[i][1]];
                tempnewnode[tnncount][1] = nodeconnections[trackIDs[i][2]][!trackIDs[i][3]];
                //adds in points from list in trackID[0] ie the first points vector 
                int tj = 0;//position in tempnewfind
                for(int j = 0; j < comboindex[trackIDs[i][0]]; j++){
                    for(int k = 0; k < *dim + 1; k++){
                        tempnewfind[tnncount][tj][k] = findpoints[trackIDs[i][0]][j][k];
                    }
                    tj++;
                }
                //next add trackID[2] ie second point vector
                for(int j = 0; j < comboindex[trackIDs[i][2]]; j++){
                    for(int k = 0; k < *dim + 1; k++){
                        tempnewfind[tnncount][tj][k] = findpoints[trackIDs[i][2]][j][k];
                    }
                    tj++;
                }
                //go onto next
                tnncount++;
            }
        }
        //next fit our new values into existing structures :)
        //aka adjust everything into the same formate but without nodes that equal two :)
        
        //first allocate for new 
        int newcombocount = *combocount - numtwo;
        int newnicount = *nicount - numtwo;
        printf("new counts => %d %d\n",newcombocount,newnicount);
        int **newnodeconnections = malloc(newcombocount * sizeof(int*));//for new nodeconnections
        int **newnodeindex = malloc(newnicount * sizeof(int*));//for new nodeindex
        int nni = 0;
        double ***newfindpoints = malloc(newcombocount * sizeof(double**));//allocs new vectors for new findpoints:)
        int nfpi = 0;//keeps track of current newfindpoints location
        int *newcomboindex = calloc(newcombocount, sizeof(int));//for new comboindex
        for(int i = 0; i < newnicount; i++){
            newnodeindex[i] = malloc(3 * sizeof(int));
        }
        for(int i = 0; i < newcombocount; i++){
            newnodeconnections[i] = calloc(2,sizeof(int));
        }
        //next go through old and add in if allowable, note we will reclassify id's :/ unfortunate but makes life easier later
        //to mitigate cringe we can add in unclassified olds and then news, as news opperate on different dimensiions anyways
        for(int i = 0; i < *combocount; i++){
            printf("loop1 i:%d, nni:%d, newnicount:%d:\n",i,nni,newnicount);
            //check both are clear to be added
            if(countIDs[nodeconnections[i][0]] != 2 && countIDs[nodeconnections[i][1]] != 2){
                //identify if need to add id0 or id1 to index/ assing existing index
                bool pass0 = true;
                int i0 = 0;
                bool pass1 = true;
                int i1 = 0;
                for(int j = 0; j < nni; j++){
                    printf("\ncomparing [%d,%d,%d] to [%d,%d,%d] \n\n",newnodeindex[j][0],newnodeindex[j][1],newnodeindex[j][2],nodeindex[nodeconnections[i][0]][0],nodeindex[nodeconnections[i][0]][1],nodeindex[nodeconnections[i][0]][2]);
                    printf("\ncomparing [%d,%d,%d] to [%d,%d,%d] \n\n",newnodeindex[j][0],newnodeindex[j][1],newnodeindex[j][2],nodeindex[nodeconnections[i][1]][0],nodeindex[nodeconnections[i][1]][1],nodeindex[nodeconnections[i][1]][2]);
                    if(newnodeindex[j][0] == nodeindex[nodeconnections[i][0]][0] && newnodeindex[j][1] == nodeindex[nodeconnections[i][0]][1] && newnodeindex[j][2] == nodeindex[nodeconnections[i][0]][2]){
                        printf("no pass0-0\n");
                        pass0 = false;//current 0 hit as already made into an id
                        i0 = j;
                    }
                    if(newnodeindex[j][0] == nodeindex[nodeconnections[i][1]][0] && newnodeindex[j][1] == nodeindex[nodeconnections[i][1]][1] && newnodeindex[j][2] == nodeindex[nodeconnections[i][1]][2]){
                        printf("no pass0-1\n");
                        pass1 = false;//current 0 hit as already made into an id
                        i1 = j;
                    }
                }
                if(pass0){
                    newnodeindex[nni][0] = nodeindex[nodeconnections[i][0]][0];
                    newnodeindex[nni][1] = nodeindex[nodeconnections[i][0]][1];
                    newnodeindex[nni][2] = nodeindex[nodeconnections[i][0]][2];
                    printf("passed0-0 with [%d,%d,%d]\n", newnodeindex[nni][0],newnodeindex[nni][1],newnodeindex[nni][2]);
                    i0 = nni;
                    nni++;
                }
                if(pass1){
                    newnodeindex[nni][0] = nodeindex[nodeconnections[i][1]][0];
                    newnodeindex[nni][1] = nodeindex[nodeconnections[i][1]][1];
                    newnodeindex[nni][2] = nodeindex[nodeconnections[i][1]][2];
                    printf("passed0-1 with [%d,%d,%d]\n", newnodeindex[nni][0],newnodeindex[nni][1],newnodeindex[nni][2]);
                    i1 = nni;
                    nni++;
                }
                //next alloc and assign find points and length, as should be 1 to 1
                newfindpoints[nfpi] = malloc(comboindex[i] * sizeof(double*));
                for(int j = 0; j < comboindex[i]; j++){
                    newfindpoints[nfpi][j] = calloc(*dim + 1,sizeof(double));
                    for(int k = 0; k < *dim + 1; k++){
                        newfindpoints[nfpi][j][k] = findpoints[i][j][k];
                    }
                }
                newcomboindex[nfpi] = comboindex[i];
                //finally set the node connections
                newnodeconnections[nfpi][0] = i0;
                newnodeconnections[nfpi][1] = i1;
                nfpi++;
            }
        }
        //next we go through and add in our temp/ combination variables 
        //we do after to apply different resctrictions
        for(int i = 0; i < numtwo; i++){
            printf("loop2 i:%d, nni:%d, newnicount:%d:\n",i,nni,newnicount);
            //tempnewnode = [numtwo][oi0,oi1]
            //we figure out the index conditions of the system and what we are adding in
            bool pass0 = true;
            int i0 = 0;
            bool pass1 = true;
            int i1 = 0;
            for(int j = 0; j < nni; j++){
                    printf("\ncomparing [%d,%d,%d] to [%d,%d,%d] \n\n",newnodeindex[j][0],newnodeindex[j][1],newnodeindex[j][2],nodeindex[tempnewnode[i][0]][0],nodeindex[tempnewnode[i][0]][1],nodeindex[tempnewnode[i][0]][2]);
                    printf("\ncomparing [%d,%d,%d] to [%d,%d,%d] \n\n",newnodeindex[j][0],newnodeindex[j][1],newnodeindex[j][2],nodeindex[tempnewnode[i][1]][0],nodeindex[tempnewnode[i][1]][1],nodeindex[tempnewnode[i][1]][2]);
                if(nodeindex[tempnewnode[i][0]][0] == newnodeindex[j][0] && nodeindex[tempnewnode[i][0]][1] == newnodeindex[j][1] && nodeindex[tempnewnode[i][0]][2] == newnodeindex[j][2]){
                    printf("no pass2-0\n");
                    pass0 = false;
                    i0 = j;
                }
                if(nodeindex[tempnewnode[i][1]][0] == newnodeindex[j][0] && nodeindex[tempnewnode[i][1]][1] == newnodeindex[j][1] && nodeindex[tempnewnode[i][1]][2] == newnodeindex[j][2]){
                    printf("no pass2-1\n");
                    pass1 = false;
                    i1 = j;
                }
            }
            if(pass0){
                newnodeindex[nni][0] = nodeindex[tempnewnode[i][0]][0];
                newnodeindex[nni][1] = nodeindex[tempnewnode[i][0]][1];
                newnodeindex[nni][2] = nodeindex[tempnewnode[i][0]][2];
                printf("passed2-0 with [%d,%d,%d]\n", newnodeindex[nni][0],newnodeindex[nni][1],newnodeindex[nni][2]);
                i0 = nni;
                nni++;
            }
            if(pass1){
                newnodeindex[nni][0] = nodeindex[tempnewnode[i][1]][0];
                newnodeindex[nni][1] = nodeindex[tempnewnode[i][1]][1];
                newnodeindex[nni][2] = nodeindex[tempnewnode[i][1]][2];
                printf("passed2-1 with [%d,%d,%d]\n", newnodeindex[nni][0],newnodeindex[nni][1],newnodeindex[nni][2]);
                i1 = nni;
                nni++;
            }

            //we can add in our points and comboindex
            newfindpoints[nfpi] = malloc(tempfindcount[i] * sizeof(double*));
            for(int j = 0; j < tempfindcount[i]; j++){
                newfindpoints[nfpi][j] = calloc(*dim + 1,sizeof(double));
                for(int k = 0; k < *dim + 1; k++){
                    newfindpoints[nfpi][j][k] = tempnewfind[i][j][k];
                }
            }
            newcomboindex[nfpi] = tempfindcount[i];
            //finally add for success
            newnodeconnections[nfpi][0] = i0;
            newnodeconnections[nfpi][1] = i1;
            nfpi++;
        }
        //now we should theoretically have everything full with the right values and sizes
        //so lastly we free current vars and setin our new ones :)
        freeCombo(&findpoints,&comboindex,&nodeconnections,&nodeindex,combocount,nicount);
        //and then we set those freed values to our new :))
        findpoints = newfindpoints;
        comboindex = newcomboindex;
        nodeconnections = newnodeconnections;
        nodeindex = newnodeindex;
        *combocount = newcombocount;
        *nicount = newnicount;
        //finally free
        //freeCombo(&newfindpoints,&newcomboindex,&newnodeconnections,&newnodeindex,&newcombocount,&newnicount);
        free(countIDs);
        for(int i = 0; i < *nicount; i++){
            free(trackIDs[i]);
        }
        free(trackIDs);
        for (int i = 0; i < numtwo; i++) {
            free(tempnewnode[i]);
            for (int j = 0; j < tempfindcount[i]; j++) {
                free(tempnewfind[i][j]);
            }
            free(tempnewfind[i]);
        }
        free(tempnewnode);
        free(tempnewfind);
        free(tempfindcount);
        //Finally we resend our structure to the function as it will check if there are any more twos needed to be merged now. 
        //This will happen until all 2's have been merged 
        
        thinNodePoint(&sDmain,dim,&findpoints,&comboindex,&nodeconnections,&nodeindex,combocount,nicount,mxpt,t);
    }//note if no numtwos then we moveon
    //push pointers
    *pfindpoints = findpoints;
    *pnodeindex = nodeindex;
    *pnodeconnections = nodeconnections;
    *pcomboindex = comboindex;
    *sd = sDmain;
}
void makeNodePoint(struct skeleDensity **sd,int *dim){
    struct skeleDensity *sDmain = *sd;
    for(int i = 0; i < *sDmain->row; i++){
        for(int j = 0; j < *sDmain->col; j++){
            for(int k = 0; k < *sDmain->dep; k++){
                if(*(sDmain->sB[i][j][k]).leng > 0){
                    *(sDmain->pcount) = *(sDmain->pcount) + 1;
                    //fprintf(stdout,"\n");
                    //fprintf(stdout,"leng1=%d\n",*(sDmain->sB[i][j][k]).leng);
                    //Here we will loop through all of our cells
                    //We need to make sure we wont scan outofbound cells
                    //So creating a check map will be beneficial
                    int **searchid = malloc(26 * sizeof(int*));//allocate 26 potential cells of reach, if a full 3x3   
                    for(int q = 0; q < 26; q++){
                        searchid[q] = calloc(3 , sizeof(int));//initalize all as [0,0,0]
                    }
                    int len = 0;
                    for(int ic = -1; ic <= 1; ic++){
                        for(int jc = -1; jc <= 1; jc++){
                            for(int kc = -1; kc <= 1; kc++){
                                int ti = i + ic;
                                int tj = j + jc;
                                int tk = k + kc;
                                if(ti >= 0 && ti < *sDmain->row){
                                    if(tj >= 0 && tj < *sDmain->col){
                                        if(tk >= 0 && tk < *sDmain->dep){
                                            //fprintf(stdout,"leng@look=%d,[%d,%d]\n",*(sDmain->sB[ti][tj][tk]).leng,ic,jc);
                                            if(!(ic == 0 && jc == 0 && kc == 0) ){// && passcheck(ic,jc,kc,dim)){
                                                if(*(sDmain->sB[ti][tj][tk]).leng != 0){
                                                    //Finally we check if the roi is connecting:)
                                                    //if(connectingROI((sDmain->sB[i][j][k]).roi,(sDmain->sB[ti][tj][tk]).roi,dim)){
                                                        searchid[len][0] = ic;  
                                                        searchid[len][1] = jc;  
                                                        searchid[len][2] = kc;
                                                        len++; 
                                                    //}
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //now we have created a list of spots we want to visit that is allowable
                    bool localmax = false;
                    bool endpoint = false;
                    //double *dense = (sDmain->sB[i][j][k]).density;
                    //for(int runi = 0; runi < len; runi++){
                    //    //goes through all our wanted cells
                    //    //fprintf(stdout,"going to look @ [%d,%d,%d]\n",i + searchid[runi][0],j + searchid[runi][1],k + searchid[runi][2]);
                    //    if(*dense < *(sDmain->sB[i+searchid[runi][0]][j+searchid[runi][1]][k+searchid[runi][2]]).density){
                    //        localmax = false;
                    //        break;
                    //    }
                    //}
                    //Here we calculate local maximum based on amount of connecting branches
                    //if(localmax){
                    //    for(int runi = 0; runi < len; runi++){
                    //        //goes through all our wanted cells
                    //        //fprintf(stdout,"going to look @ [%d,%d,%d]\n",i + searchid[runi][0],j + searchid[runi][1],k + searchid[runi][2]);
                    //        if(len < *(sDmain->sB[i+searchid[runi][0]][j+searchid[runi][1]][k+searchid[runi][2]]).connections){
                    //            localmax = false;
                    //            break;
                    //        }
                    //    }
                    //}
                    //if(len <= 2){
                    //    localmax = false;
                    //}
                    if(len == 1  || len == 0){
                        endpoint = true;
                    }
                    else if(*dim == 2 && len == 2){
                        //Here we check if the two are touching eachother 
                        //we are also in 2D so it will be simpler
                        int r1 = searchid[0][0];
                        int r2 = searchid[1][0];
                        int c1 = searchid[0][1];
                        int c2 = searchid[1][1];
                        if((r1 == r2 && (abs(c1-c2) == 1)) || (c1 == c2 && (abs(r1-r2) == 1))){
                            endpoint = true;
                        }
                    }
                    else if(*dim == 2 && len == 3){
                        //Here we check if there is a corner, and the two neighboring cells
                        int r1 = searchid[0][0];
                        int r2 = searchid[1][0];
                        int r3 = searchid[2][0];
                        int c1 = searchid[0][1];
                        int c2 = searchid[1][1];
                        int c3 = searchid[2][1];
                        //first check if one node is a corner
                        int ccase = 0;
                        if(r1 != 0 && c1 != 0){
                            ccase += 1;
                        }
                        if(r2 != 0 && c2 != 0){
                            ccase += 10;
                        }
                        if(r3 != 0 && c3 != 0){
                            ccase += 100;
                        }
                        //Now we have a index of the corner
                        //We need to check if the other two are each touching in opposite dimensions
                        int nr1;
                        int nr2;
                        int nc1;
                        int nc2;
                        if(ccase == 1){
                            nr1 = searchid[1][0];
                            nr2 = searchid[2][0];
                            nc1 = searchid[1][1];
                            nc2 = searchid[2][1];
                            ccase = 0;
                            }
                        else if(ccase == 10){
                            nr1 = searchid[0][0];
                            nr2 = searchid[2][0];
                            nc1 = searchid[0][1];
                            nc2 = searchid[2][1];
                            ccase = 1;
                        }
                        else if(ccase == 100){
                            nr1 = searchid[0][0];
                            nr2 = searchid[1][0];
                            nc1 = searchid[0][1];
                            nc2 = searchid[1][1];
                            ccase = 2;
                        }
                        else{
                            ccase = -1;
                        }
                        if(ccase != -1){
                            if((searchid[ccase][0] == nr1 && searchid[ccase][1] == nc2) || (searchid[ccase][0] == nr2 && searchid[ccase][1] == nc1)){
                                endpoint = true;
                            }
                        }
                        else{
                            if((r1 == r2 && r1 == r3) || (c1 == c2 && c1 == c3)){
                                endpoint = true;
                            }
                        }
                    }
                    if(!endpoint && *dim == 2 && len > 0){
                        int corners = 0;
                        int sides = 0;
                        bool xmax = false;
                        bool ymax = false;
                        bool xmin = false;
                        bool ymin = false;
                        for(int runi = 0; runi < len; runi++){
                            if(searchid[runi][0] != 0 && searchid[runi][1] != 0){
                                corners++;
                            }
                            if((searchid[runi][0] != 0 && searchid[runi][1] == 0) || (searchid[runi][0] == 0 && searchid[runi][1] != 0)){
                                sides++;
                            }
                            if(searchid[runi][0] == 1){
                                xmax = true;
                            }
                            else if(searchid[runi][0] == -1){
                                xmin = true;
                            }
                            if(searchid[runi][1] == 1){
                                ymax = true;
                            }
                            else if(searchid[runi][1] == -1){
                                ymin = true;
                            }
                        }
                        if(corners+sides >= 4 && (xmax && ymax && xmin && ymin)){
                            if(corners+sides == 4){
                                //Checks that the 4 isnt crossing through everything 
                                if(corners > 2 || sides > 2){
                                    localmax = true;
                                }
                                else{
                                    //sometimes there is a case where sides = 2 &  corners = 2 & it is a bifurcation area
                                    //so we target that configuration here 
                                    int *corner1 = calloc(2,sizeof(int));
                                    int *corner2 = calloc(2,sizeof(int));
                                    int *side1 = calloc(2,sizeof(int));
                                    int *side2 = calloc(2,sizeof(int));
                                    bool addc = false;
                                    bool adds = false;
                                    for(int runi = 0; runi < len; runi++){
                                        //fprintf(stdout,"searchid[%d,%d]\n",searchid[runi][0],searchid[runi][1]);
                                        if(searchid[runi][0] != 0 && searchid[runi][1] != 0){
                                            if(!addc){
                                                corner1[0] = searchid[runi][0];
                                                corner1[1] = searchid[runi][1];
                                                addc = true;
                                            }
                                            else{
                                                corner2[0] = searchid[runi][0];
                                                corner2[1] = searchid[runi][1];
                                            }
                                        }
                                        if((searchid[runi][0] != 0 && searchid[runi][1] == 0) || (searchid[runi][0] == 0 && searchid[runi][1] != 0)){
                                            if(!adds){
                                                side1[0] = searchid[runi][0];
                                                side1[1] = searchid[runi][1];
                                                adds = true;
                                            }
                                            else{
                                                side2[0] = searchid[runi][0];
                                                side2[1] = searchid[runi][1];
                                            }
                                        }
                                    }
                                    //fprintf(stdout,"try correction [%d,%d],[%d,%d]\n",corner1[0],corner1[1],corner2[0],corner2[1]);
                                    if(corner1[0] == corner2[0] || corner1[1] == corner2[1]){
                                        //fprintf(stdout,"correction\n");
                                        localmax = true;
                                    }
                                    free(corner1);
                                    free(corner2);
                                    free(side1);
                                    free(side2);
                                }
                            }
                            else{
                                localmax = true;
                            }
                        }
                    }
                    if(endpoint){
                        //If it is a end point, we will make a node point of the furthest point :)
                        //fprintf(stdout,"endpoint\n");
                        //fprintf(stdout,"@[%d,%d]\n",i,j);
                        //fprintf(stdout,"len = %d\n",len);
                        //for(int q = 0; q < len; q++){
                        //    fprintf(stdout,"connections = [%d,%d]\n",searchid[q][0],searchid[q][1]);
                        //}
                        struct skeleBounds *localcell = &(sDmain->sB[i][j][k]);
                        *localcell->hasNode = true;
                        double extra = 1;
                        double *calcpoint =  calloc(*dim+extra,sizeof(double));
                        //Now we have our space malloced
                        //Next we calculate our Node Point based on the average of the cell
                        //To find our wanted node, we take the center of our data, which is calulated in the skeleBounds
                        double cx = (*(sDmain->xmin) + *(sDmain->xmax)) / 2;
                        double cy = (*(sDmain->ymin) + *(sDmain->ymax)) / 2;
#if dimension == 3
                        double cz = (*(sDmain->zmin) + *(sDmain->zmax)) / 2;
#endif
                        //And we will select the furthest point from the 'center'
                        double furthest = 0.;
                        int holdq = 0;
                        for(int q = 0; q < *(localcell->leng); q++){
                            double dx = localcell->points[q][0] - cx;
                            double dy = localcell->points[q][1] - cy;
#if dimension == 2
                            double d = pow(dx,2) + pow(dy,2); 
#else
                            double dz = localcell->points[q][2] - cz;
                            double d = pow(dx,2) + pow(dy,2) + pow(dz,2); 
#endif
                            d = sqrt(d);
                            if(d > furthest){
                                furthest = d;
                                holdq = q;
                            }
                        }
                        calcpoint[0] = localcell->points[holdq][0];
                        calcpoint[1] = localcell->points[holdq][1];
                        if(*dim == 3){
                            calcpoint[2] = localcell->points[holdq][2];
                            calcpoint[3] = localcell->points[holdq][3];
                        }
                        else{
                            calcpoint[2] = localcell->points[holdq][2];
                        }
                        localcell->nodepoint = calcpoint;
                        *localcell->smode = 2;
                        (*(sDmain->ncount))++;
                    }
                    else if(localmax){
                        //fprintf(stdout,"localmax\n");
                        //fprintf(stdout,"@[%d,%d]\n",i,j);
                        //fprintf(stdout,"len = %d\n",len);
                        //for(int q = 0; q < len; q++){
                        //    fprintf(stdout,"connections = [%d,%d]\n",searchid[q][0],searchid[q][1]);
                        //}
                        //If it is a local maximum, we will make a node point :)
                        struct skeleBounds *localcell = &(sDmain->sB[i][j][k]);
                        *localcell->hasNode = true;
                        double extra = 1;
                        double *calcpoint =  calloc(*dim+extra,sizeof(double));
                        //Now we have our space malloced
                        //Next we calculate our Node Point based on the average of the cell
                        double ax = 0.;
                        double ay = 0.;
                        double ar = 0.;
#if dimension == 3
                        double az = 0.;
#endif
                        //fprintf(stdout,"leng=%d\n",*(localcell->leng));
                        for(int q = 0; q < *(localcell->leng); q++){
                            ax += localcell->points[q][0];
                            ay += localcell->points[q][1];
#if dimension == 2
                            ar += localcell->points[q][2];
#else
                            az += localcell->points[q][2];
                            ar += localcell->points[q][3];
#endif
                            //fprintf(stdout,"summing[%f,%f,%f]\n",ax,ay,ar);
                        }
                        ax = ax / *(localcell->leng);
                        ay = ay / *(localcell->leng);
#if dimension == 3
                        az = az / *(localcell->leng);
#endif      
                        ar = ar / *(localcell->leng);
                        calcpoint[0] = ax;
                        calcpoint[1] = ay;
#if dimension == 3    
                        calcpoint[2] = az;
                        calcpoint[3] = ar;
#else
                        calcpoint[2] = ar;
#endif      
                        localcell->nodepoint = calcpoint;
                        *localcell->smode = 1;
                        (*(sDmain->ncount))++;
                    }
                    //finally free up our searching
                    for(int q = 0; q < 26; q++){
                        free(searchid[q]);
                    }
                    free(searchid);
                }
            }
        }
    }
}
//method for tracing along our skeleDensity struct 
int floodDensity(struct skeleDensity **sD,int *dim,int *nodeLocation,int distance,int ***pvisitStack,int *visitcount,int **markerids,int *lengmarker,double *tolerance){
    //Firstly we set our curent cell's distance & assign it a new 
    struct skeleDensity *sd = *sD;
    int **visitStack = *pvisitStack; 
    if(distance == 0){
        visitStack[0][0] = nodeLocation[0];
        visitStack[0][1] = nodeLocation[1];
        visitStack[0][2] = nodeLocation[2];
        for(int i = 0; i < *lengmarker; i++){
            visitStack[i+1][0] = markerids[i][0];
            visitStack[i+1][1] = markerids[i][1];
            visitStack[i+1][2] = markerids[i][2];
        }
    }
    int i = visitStack[*visitcount][0];
    int j = visitStack[*visitcount][1];
    int k = visitStack[*visitcount][2];
    //fprintf(stdout,"@[%d,%d,%d];distance=%d,vc=%d\n",i,j,k,distance,*visitcount);
    if(*(sd->sB[i][j][k]).closedis1 == -1 || *(sd->sB[i][j][k]).closedis1 > distance){
        //*(sd->sB[i][j][k]).closedis2 = *(sd->sB[i][j][k]).closedis1;
        //(sd->sB[i][j][k]).closeid2[0] = (sd->sB[i][j][k]).closeid1[0];
        //(sd->sB[i][j][k]).closeid2[1] = (sd->sB[i][j][k]).closeid1[1];
        //(sd->sB[i][j][k]).closeid2[2] = (sd->sB[i][j][k]).closeid1[2];
        *(sd->sB[i][j][k]).closedis1 = distance;
        (sd->sB[i][j][k]).closeid1[0] = nodeLocation[0];
        (sd->sB[i][j][k]).closeid1[1] = nodeLocation[1];
        (sd->sB[i][j][k]).closeid1[2] = nodeLocation[2];
        *visitcount = *visitcount + +1;
        //fprintf(stdout,"closest node[%d,%d,%d]\n",(sd->sB[i][j][k]).closeid1[0],(sd->sB[i][j][k]).closeid1[1],(sd->sB[i][j][k]).closeid1[2]);
    }
    else{
        //fprintf(stdout,"pass on node, belongs to[%d,%d,%d] (%d vs %d)\n",(sd->sB[i][j][k]).closeid1[0],(sd->sB[i][j][k]).closeid1[1],(sd->sB[i][j][k]).closeid1[2],distance,*(sd->sB[i][j][k].closedis1));
        //if(*(sd->sB[i][j][k]).closedis2 == -1 || *(sd->sB[i][j][k]).closedis2 > distance){
        //    *(sd->sB[i][j][k]).closedis2 = distance;
        //    (sd->sB[i][j][k]).closeid2[0] = nodeLocation[0];
        //    (sd->sB[i][j][k]).closeid2[1] = nodeLocation[1];
        //    (sd->sB[i][j][k]).closeid2[2] = nodeLocation[2];
        //}
        //else{
            return 0;
        //}
    }
    //Next we determine if there are any moves we want to make
    int **searchid = malloc(26 * sizeof(int*));//allocate 26 potential cells of reach, if a full 3x3   
    for(int q = 0; q < 26; q++){
        searchid[q] = calloc(3 , sizeof(int));//initalize all as [0,0,0]
    }
    int len = 0;
    for(int ic = -1; ic <= 1; ic++){
        for(int jc = -1; jc <= 1; jc++){
            for(int kc = -1; kc <= 1; kc++){
                int ti = i + ic;
                int tj = j + jc;
                int tk = k + kc;
                if(ti >= 0 && ti < *sd->row){
                    if(tj >= 0 && tj < *sd->col){
                        if(tk >= 0 && tk < *sd->dep){
                            if(!(ic == 0 && jc == 0 && kc == 0)){
                                if(*(sd->sB[ti][tj][tk]).leng != 0){
                                    //Finally we check if the roi is connecting:)
                                    //if(connectingROI((sd->sB[i][j][k]).roi,(sd->sB[ti][tj][tk]).roi,dim,tolerance)){
                                        searchid[len][0] = ic;  
                                        searchid[len][1] = jc;  
                                        searchid[len][2] = kc;
                                        len++; 
                                    //}
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //Now we have all our search id's, we first send floods towards the directly touching pieces, unless theres non in whichcase we will 
    int edgesvisited = 0;//we track visited edges, and in cases where all have been gone to, we will bridge to a corner
    int corners = 0;//we track corners/edges/cases where more than 1 variable isnt 0,
    int edges = 0;
    //this will let us know if we need to reach to corners or not
    //fprintf(stdout,"\n");
    for(int q = 0; q < len; q++){
        //fprintf(stdout,"trying for [%d,%d,%d] (%d/%d)\n",i+searchid[q][0],j+searchid[q][1],k+searchid[q][2],q+1,len);
        if((searchid[q][0] != 0 && searchid[q][1] == 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] != 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] == 0 && searchid[q][2] != 0)){
            //bool allowpass = true;
            //for(int p = 0; p < *visitcount; p++){
            //    //check if node is already visited
            //    //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
            //    if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
            //        allowpass = false;
            //        break;
            //    }
            //}
            //if(allowpass){
                //Here we have determined that we are next to the cell & it hasnt been visited yet, so we will send a flood
                //fprintf(stdout,"going to%d++\n",*visitcount);
                visitStack[*visitcount][0] = i + searchid[q][0];
                visitStack[*visitcount][1] = j + searchid[q][1];
                visitStack[*visitcount][2] = k + searchid[q][2];
                //*visitcount = *visitcount + 1;
                edgesvisited = edgesvisited + floodDensity(&sd,dim,nodeLocation,distance + 1,&visitStack,visitcount,markerids,lengmarker,tolerance);
                edges++;
            //}
        }
        else{
            //fprintf(stdout,"wascorner[%d,%d,%d]\n",searchid[q][0],searchid[q][1],searchid[q][2]);
            //bool allowpass = true;
            //for(int p = 0; p < *visitcount; p++){
            //    //check if node is already visited
            //    //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
            //    if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
            //        allowpass = false;
            //        break;
            //    }
            //}
            //if(allowpass){
                corners++;
            //}
        }
    }
    //now we correct to corners if no other options
    if(corners > 0){
        for(int q = 0; q < len; q++){
            //fprintf(stdout,"trying for [%d,%d,%d] (%d/%d)\n",i+searchid[q][0],j+searchid[q][1],k+searchid[q][2],q+1,len);
            if((searchid[q][0] != 0 && searchid[q][1] != 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] != 0 && searchid[q][2] != 0) || (searchid[q][0] != 0 && searchid[q][1] == 0 && searchid[q][2] != 0)){
                //bool allowpass = true;
                //for(int p = 0; p < *visitcount; p++){
                //    //check if node is already visited
                //    //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
                //    if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
                //        allowpass = false;
                //        break;
                //    }
                //}
                //if(allowpass){
                    //Here we have determined that we are next to the cell & it hasnt been visited yet, so we will send a flood
                    //fprintf(stdout,"going to%d++\n",*visitcount);
                    visitStack[*visitcount][0] = i + searchid[q][0];
                    visitStack[*visitcount][1] = j + searchid[q][1];
                    visitStack[*visitcount][2] = k + searchid[q][2];
                    //*visitcount = *visitcount + 1;
                    if(edges == 1){
                        floodDensity(&sd,dim,nodeLocation,distance + 1,&visitStack,visitcount,markerids,lengmarker,tolerance);
                    }
                    else{
                        floodDensity(&sd,dim,nodeLocation,distance + 2,&visitStack,visitcount,markerids,lengmarker,tolerance);
                    }
                //}
            }
        }
    }
    //fprintf(stdout,"\n");
    for(int q = 0; q < 26; q++){
        free(searchid[q]);
    }
    free(searchid);
    *pvisitStack = visitStack;
    *sD = sd;
    return 1;
}
//resets all visited stacks to 0
void resetFloodStack(int **vs, int leng, int **ids, int idleng){
    for(int i = 0; i < leng; i++){
        vs[i][0] = 0;
        vs[i][1] = 0;
        vs[i][2] = 0;
    }
    for(int i = 0; i < idleng; i++){
        ids[i][0] = 0;
        ids[i][1] = 0;
        ids[i][2] = 0;
    }

}
//Method for reducing LocalMax's where more than one exist in a 3x3 around it, uses endpoint calcs to 
//optimize placement
void reduceLocalMax(struct skeleDensity **sD,int *dim){
    struct skeleDensity *sd  = *sD;
    for(int i = 0; i < *sd->row; i++){
        for(int j = 0; j < *sd->col; j++){
            for(int k = 0; k < *sd->dep; k++){
                if(*(sd->sB[i][j][k]).hasNode && *(sd->sB[i][j][k]).smode  == 1){
                    //We have reached a point which is considered a localmax point 
                    //and now we check if there are any point around it
                    int tcount = 0;
                    int tendcount = 0;
                    bool closecase = false;
                    //works for a 5x5 grid
                    for(int ic = -2; ic <= 2; ic++){
                        int ti = i + ic;
                        for(int jc = -2; jc <= 2; jc++){
                            int tj = j + jc;
                            for(int kc = -2; kc <= 2; kc++){
                                int tk = k + kc;
                                if(ti >= 0 && ti < *sd->row){
                                    if(tj >= 0 && tj < *sd->col){
                                        if(tk >= 0 && tk < *sd->dep){
                                            if(!(ic == 0 && jc == 0 && kc == 0)){
                                                //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                //Next we shall count the amount of nearby nodes made in smode = 1
                                                if(*(sd->sB[ti][tj][tk]).smode == 2){
                                                    closecase = true;
                                                    tendcount++;
                                                }
                                                if(*(sd->sB[ti][tj][tk]).smode == 1){
                                                    tcount++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //Now we have count & if its too close to a endpoint.
                    if(closecase){
                        for(int ic = -2; ic <= 2; ic++){
                            int ti = i + ic;
                            for(int jc = -2; jc <= 2; jc++){
                                int tj = j + jc;
                                for(int kc = -2; kc <= 2; kc++){
                                    int tk = k + kc;
                                    if(ti >= 0 && ti < *sd->row){
                                        if(tj >= 0 && tj < *sd->col){
                                            if(tk >= 0 && tk < *sd->dep){
                                                //if(!(ic == 0 && jc == 0 && kc == 0)){
                                                    //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                    //Next we shall count the amount of nearby nodes made in smode = 1
                                                    if(*(sd->sB[ti][tj][tk]).smode == 1){
                                                        //because too close to endpoint, we disable the nodepoint
                                                        *(sd->sB[ti][tj][tk]).hasNode = false;
                                                        free((sd->sB[ti][tj][tk]).nodepoint);
                                                        *(sd->sB[ti][tj][tk]).smode = 0;
                                                        *(sd->ncount) = *(sd->ncount) - 1;
                                                    }
                                                //}
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //Lastly we will also merge close end points, these should be okay
                        if(tendcount > 1){
                            int hcount = 0;
                            for(int ic = -2; ic <= 2; ic++){
                                for(int jc = -2; jc <= 2; jc++){
                                    for(int kc = -2; kc <= 2; kc++){
                                        int ti = i + ic;
                                        int tj = j + jc;
                                        int tk = k + kc;
                                        if(ti >= 0 && ti < *sd->row){
                                            if(tj >= 0 && tj < *sd->col){
                                                if(tk >= 0 && tk < *sd->dep){
                                                    if(!(ic == 0 && jc == 0 && kc == 0)){
                                                        //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                        //Next we shall count the amount of nearby nodes made in smode = 1
                                                        if(*(sd->sB[ti][tj][tk]).smode == 2){
                                                            //because too close to nodepoint, we disable the nodepoint
                                                            if(hcount == 0){
                                                                hcount++;
                                                            }
                                                            else{
                                                                *(sd->sB[ti][tj][tk]).hasNode = false;
                                                                free((sd->sB[ti][tj][tk]).nodepoint);
                                                                *(sd->sB[ti][tj][tk]).smode = 0;
                                                                *(sd->ncount) = *(sd->ncount) - 1;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //And if not end point, it will collect a new nodepoint at the node nearby with the mode around it
                    else if(tcount > 0){
                        //For choosing which nodepoint to save, we will first average them up, and then 
                        //given it falls inside one of the bounds of a point, that will be considered the branch
                        int iavg = 0;
                        int javg = 0;
                        int kavg = 0;
                        int avgcnt = 0;
                        for(int ic = -2; ic <= 2; ic++){
                            int ti = i + ic;
                            for(int jc = -2; jc <= 2; jc++){
                                int tj = j + jc;
                                for(int kc = -2; kc <= 2; kc++){
                                    int tk = k + kc;
                                    if(ti >= 0 && ti < *sd->row){
                                        if(tj >= 0 && tj < *sd->col){
                                            if(tk >= 0 && tk < *sd->dep){
                                                //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                //Next we shall count the amount of nearby nodes made in smode = 1
                                                if(*(sd->sB[ti][tj][tk]).smode == 1){
                                                    iavg = iavg + ti;
                                                    javg = javg + tj;
                                                    kavg = kavg + tk;
                                                    avgcnt++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        //Now we average the amount
                        iavg = (int)round(iavg / avgcnt);
                        javg = (int)round(javg / avgcnt);
                        kavg = (int)round(kavg / avgcnt);
                        for(int ic = -2; ic <= 2; ic++){
                            int ti = i + ic;
                            for(int jc = -2; jc <= 2; jc++){
                                int tj = j + jc;
                                for(int kc = -2; kc <= 2; kc++){
                                    int tk = k + kc;
                                    if(ti >= 0 && ti < *sd->row){
                                        if(tj >= 0 && tj < *sd->col){
                                            if(tk >= 0 && tk < *sd->dep){
                                                //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                //Next we shall count the amount of nearby nodes made in smode = 1
                                                if(*(sd->sB[ti][tj][tk]).smode == 1 && !(ti == iavg && tj == javg && tk == kavg)){
                                                    *(sd->sB[ti][tj][tk]).hasNode = false;
                                                    free((sd->sB[ti][tj][tk]).nodepoint);
                                                    *(sd->sB[ti][tj][tk]).smode = 0;
                                                    *(sd->ncount) = *(sd->ncount) - 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }   
    }
    *sD = sd;
}
//handles allocation/reallocation of combos
void addCombo(int *masterid,int *indexCombo,int *comboCount,int ***pnodeconnections,int ***pnodeindex,int *nicount){
    int **nodeconnections = *pnodeconnections;
    int **nodeindex = *pnodeindex;
    int im = -1;
    int ic = -1;
    //find index of id's if already exist
    for(int i = 0; i < *nicount; i++){
        if((masterid[0]==nodeindex[i][0]) && (masterid[1]==nodeindex[i][1]) && (masterid[2]==nodeindex[i][2])){
            im = i;
        }
        else if((indexCombo[0]==nodeindex[i][0]) && (indexCombo[1]==nodeindex[i][1]) && (indexCombo[2]==nodeindex[i][2])){
            ic = i;
        }
    }
    //First calc needed alloc space
    int nialloc = 0;
    if(im == -1){
        printf("im\n");
        nialloc++;
    }
    if(ic == -1){
        printf("ic\n");
        nialloc++;
    }
    int **tempnodeindex = malloc((nialloc + *nicount) * sizeof(int*));
    int mode = 0;
    //create new nodeindex
    for(int i = 0; i < (nialloc + *nicount); i++){
        tempnodeindex[i] = calloc(3,sizeof(int));
        if(i < *nicount){
            tempnodeindex[i][0] = nodeindex[i][0];
            tempnodeindex[i][1] = nodeindex[i][1];
            tempnodeindex[i][2] = nodeindex[i][2];
        }
        else{
            if(im == -1 && mode == 0){
                mode++;
                tempnodeindex[i][0] = masterid[0];
                tempnodeindex[i][1] = masterid[1];
                tempnodeindex[i][2] = masterid[2];
                im = i;
            }
            else if(ic == -1){
                tempnodeindex[i][0] = indexCombo[0];
                tempnodeindex[i][1] = indexCombo[1];
                tempnodeindex[i][2] = indexCombo[2];
                ic = i;
            }
        }
    }
    //free old
    for(int i = 0; i < *nicount; i++){
        free(nodeindex[i]);
    }
    if(nodeindex != NULL){
        free(nodeindex);
    }
    //set new
    *pnodeindex = tempnodeindex;
    *nicount = *nicount + nialloc;
    if(im == -1 || ic == -1){
        printf("comboadd error !\n");
    }
    //Next we add the combo with the respective index
    int **tempnodeconnections = malloc((*comboCount + 1) * sizeof(int*));
    for(int i = 0; i < (*comboCount + 1); i++){
        tempnodeconnections[i] = calloc(2,sizeof(int));
        if(i < *comboCount){
            tempnodeconnections[i][0] = nodeconnections[i][0];
            tempnodeconnections[i][1] = nodeconnections[i][1];
        }
        else{
            tempnodeconnections[i][0] = im;
            tempnodeconnections[i][1] = ic;
        }
    }
    for(int i = 0; i < *comboCount; i++){
        free(nodeconnections[i]);
    }
    if(nodeconnections != NULL){
        free(nodeconnections);
    }
    *pnodeconnections = tempnodeconnections;
    *comboCount = *comboCount + 1;
}
//floodpart2 for checking id mesh
void floodDensity2(struct skeleDensity **sD,int *dim,int *combocount,int ***pnodeconnections,int ***pnodeindex,int *nicount,int **pmasterid,int ***pgonestack,int *gonecount){
    struct skeleDensity *sd  = *sD;
    int **nodeconnections = *pnodeconnections;
    int **nodeindex = *pnodeindex;
    int **gonestack = *pgonestack;
    int *masterid = *pmasterid;
    //get curent location
    int i = gonestack[*gonecount][0];
    int j = gonestack[*gonecount][1];
    int k = gonestack[*gonecount][2];
    if((sd->sB[i][j][k]).closeid1[0] == masterid[0] && (sd->sB[i][j][k]).closeid1[1] == masterid[1] && (sd->sB[i][j][k]).closeid1[2] == masterid[2]){
        //If our id is the same, we move on :)
        //first test for edge nodes & move
        for(int q = -1; q <= 1; q++){
            for(int p = -1; p <= 1; p++){
                for(int r = -1; r <= 1; r++){
                    if(!(q == 0 && p == 0 && r == 0)){
                        int ic = i + q;
                        int jc = j + p;
                        int kc = k + r;
                        if(ic >= 0 && ic < *(sd->row)){
                            if(jc >= 0 && jc < *(sd->col)){
                                if(kc >= 0 && kc < *(sd->dep)){
                                    if(*(sd->sB[ic][jc][kc]).leng != 0){
                                        //Within bounds
                                        if((q != 0 && p == 0 && r == 0) || (q == 0 && p != 0 && r == 0) || (q == 0 && p == 0 && r != 0)){
                                            //Edge node we visit first
                                            bool passcase = true;
                                            for(int l = 0; l < *gonecount; l++){
                                                if((ic == gonestack[l][0]) && (jc == gonestack[l][1]) && (kc == gonestack[l][2])){
                                                    passcase = false;
                                                    break;
                                                }
                                            }
                                            if(passcase){
                                                *gonecount = *gonecount + 1;
                                                gonestack[*gonecount][0] = ic; 
                                                gonestack[*gonecount][1] = jc; 
                                                gonestack[*gonecount][2] = kc; 
                                                floodDensity2(sD,dim,combocount,&nodeconnections,&nodeindex,nicount,&masterid,&gonestack,gonecount);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //Next test for corners and move if not visited
        for(int q = -1; q <= 1; q++){
            for(int p = -1; p <= 1; p++){
                for(int r = -1; r <= 1; r++){
                    if(!(q == 0 && p == 0 && r == 0)){
                        int ic = i + q;
                        int jc = j + p;
                        int kc = k + r;
                        if(ic >= 0 && ic < *(sd->row)){
                            if(jc >= 0 && jc < *(sd->col)){
                                if(kc >= 0 && kc < *(sd->dep)){
                                    if(*(sd->sB[ic][jc][kc]).leng != 0){
                                        //Within bounds
                                        if((q != 0 && p != 0 && r != 0) || (q != 0 && p != 0 && r == 0) || (q != 0 && p == 0 && r != 0) || (q == 0 && p != 0 && r != 0)){
                                            //now corners 
                                            bool passcase = true;
                                            for(int l = 0; l < *gonecount; l++){
                                                if((ic == gonestack[l][0]) && (jc == gonestack[l][1]) && (kc == gonestack[l][2])){
                                                    passcase = false;
                                                    break;
                                                }
                                            }
                                            if(passcase){
                                                *gonecount = *gonecount + 1;
                                                gonestack[*gonecount][0] = ic; 
                                                gonestack[*gonecount][1] = jc; 
                                                gonestack[*gonecount][2] = kc; 
                                                floodDensity2(sD,dim,combocount,&nodeconnections,&nodeindex,nicount,&masterid,&gonestack,gonecount);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else{// if(!((sd->sB[i][j][k]).closeid1[0] == 0 && (sd->sB[i][j][k]).closeid1[1] == 0 && (sd->sB[i][j][k]).closeid1[2] == 0)){
        //If its different than we check the combo count
        int *indexCombo = malloc(3*sizeof(int));
        indexCombo[0] = (sd->sB[i][j][k]).closeid1[0];
        indexCombo[1] = (sd->sB[i][j][k]).closeid1[1];
        indexCombo[2] = (sd->sB[i][j][k]).closeid1[2];
        bool clearindexcombo = false;
        //fprintf(stdout,"--combocount--%d\n",*combocount);
        for(int q = 0; q < *combocount; q++){
            if((masterid[0]==nodeindex[nodeconnections[q][0]][0]) && (masterid[1]==nodeindex[nodeconnections[q][0]][1]) && (masterid[2]==nodeindex[nodeconnections[q][0]][2])){
                if((indexCombo[0]==nodeindex[nodeconnections[q][1]][0]) && (indexCombo[1]==nodeindex[nodeconnections[q][1]][1]) && (indexCombo[2]==nodeindex[nodeconnections[q][1]][2])){
                    //Here we have found one variant
                    clearindexcombo = true;
                }

            }
            else if((masterid[0]==nodeindex[nodeconnections[q][1]][0]) && (masterid[1]==nodeindex[nodeconnections[q][1]][1]) && (masterid[2]==nodeindex[nodeconnections[q][1]][2])){
                if((indexCombo[0]==nodeindex[nodeconnections[q][0]][0]) && (indexCombo[1]==nodeindex[nodeconnections[q][0]][1]) && (indexCombo[2]==nodeindex[nodeconnections[q][0]][2])){
                    clearindexcombo = true;
                }
            }
        }
        if(clearindexcombo){
            free(indexCombo);
        }
        else{
            //We didnt find the index combo, so we will allocate for a new index combo, which places the masterid in position 1
            addCombo(masterid,indexCombo,combocount,&nodeconnections,&nodeindex,nicount);
            free(indexCombo);
        }
    }
    *pgonestack = gonestack;
    *pmasterid = masterid;
    *pnodeconnections = nodeconnections;
    *pnodeindex = nodeindex;
    *sD = sd;
}
void sortAng(double **pangles,int ***paltid,int **paltloc,int *leng,double t){
    double *angles = *pangles;
    int **altid = *paltid;
    int *altloc = *paltloc;
    bool didswap;
    for(int i = 0; i < *leng - 1;i++){
        didswap = false;
        for(int j = 0; j < *leng - i - 1;j++){
            if(angles[j] > angles[j + 1]){
                //swaps order
                //fprintf(stdout,"b4[%f-%f] , [%d,%d,%d] - [%d,%d,%d]\n",angles[j],angles[j + 1],altid[j][0],altid[j][1],altid[j][2],altid[j+1][0],altid[j+1][1],altid[j+1][2]);
                double tempa = angles[j];
                angles[j] = angles[j + 1];
                angles[j + 1] = tempa;
                int *tempid = altid[j];
                altid[j] = altid[j + 1];
                altid[j + 1] = tempid;
                int temploc = altloc[j];
                altloc[j] = altloc[j + 1];
                altloc[j + 1] = temploc;
                didswap = true;
            }
        }
        if(!didswap){
            break;
        }
    }
    *pangles = angles;
    *paltid = altid;
    *paltloc = altloc;
}
//Pathing Secondary Nodes
void sortConnections(struct skeleDensity **sD,int *dim,int *combocount,int ***pnodeconnections,double ****pfindpoints,int **pcomboindex,int ***pnodeindex,int *nicount,double t){
    struct skeleDensity *sd = *sD;
    double ***findpoints = *pfindpoints;
    int **nodeindex = *pnodeindex;
    int **nodeconnections = *pnodeconnections;
    int *comboindex = *pcomboindex;
    //First we go through and count up all the occurances of each id
    int *idcount = calloc(*nicount , sizeof(int));
    int *sidcount = calloc(*nicount , sizeof(int));
    for(int i = 0; i < *combocount; i++){
       idcount[nodeconnections[i][0]]++;//idcount counts times each id appears
       idcount[nodeconnections[i][1]]++;
       sidcount[nodeconnections[i][0]] = i;//note sidcount only helps in situations with each spline containing one idcount
       sidcount[nodeconnections[i][1]] = i;
    }
    for(int p = 0; p < *nicount; p++){
        int cid = idcount[p];
        if(cid == 1){
            //Only one occurance, so we will add points & 
            for(int i = 0 ; i < *sd->row; i++){
                for(int j = 0 ; j < *sd->col; j++){
                    for(int k = 0 ; k < *sd->dep; k++){
                        //goes through each cell && if == index, then we pass & add :)
                        if(*(sd->sB[i][j][k].leng) != 0){
                            if(((sd->sB[i][j][k]).closeid1[0]==nodeindex[p][0])&&((sd->sB[i][j][k]).closeid1[1]==nodeindex[p][1])&&((sd->sB[i][j][k]).closeid1[2] == nodeindex[p][2])){
                                //Here we now are in a cell which has the same node id as our one count, so we go through each point & assign it :)
                                for(int q = 0; q < *(sd->sB[i][j][k]).leng; q++){
                                    findpoints[sidcount[p]][comboindex[sidcount[p]]] = (sd->sB[i][j][k]).points[q];
                                    //fprintf(stdout,"added #%d to [%d,%d,%d] - [%f,%f|%f]\n",comboindex[sidcount[p]],nodeindex[p][0],nodeindex[p][1],nodeindex[p][2],findpoints[sidcount[p]][comboindex[sidcount[p]]][0],findpoints[sidcount[p]][comboindex[sidcount[p]]][1],findpoints[sidcount[p]][comboindex[sidcount[p]]][2]);
                                    comboindex[sidcount[p]]++;
                                }
                            }
                        }
                    }
                }
            }
        }
        else{
            //Multiple occurances, so we will use angles to sort id's
            //First we roughly calculate the angle based on node points
            int i = nodeindex[p][0];
            int j = nodeindex[p][1];
            int k = nodeindex[p][2];
            int ti = 0;
            int tj = 0;
            int tk = 0;
            double *curpt = (sd->sB[i][j][k]).nodepoint;//current point location
            double *comppt;//comparison point of locations
            double *angles = malloc(cid * sizeof(double));
            int **altid = malloc(cid * sizeof(int*));
            int *altloc = malloc(cid * sizeof(int));
            for(int q = 0; q < cid;q++){
                altid[q] = calloc(3,sizeof(int));
                //Since there are 'cid' connections, we calculate that many angles
                int helpcount = 0;
                for(int itemp = 0; itemp < *combocount; itemp++){
                    if(nodeconnections[itemp][0] == p){
                        if(helpcount == q){
                            altid[q][0] = nodeindex[nodeconnections[itemp][1]][0];
                            altid[q][1] = nodeindex[nodeconnections[itemp][1]][1];
                            altid[q][2] = nodeindex[nodeconnections[itemp][1]][2];
                            altloc[q] = itemp;
                            break;
                        }
                        helpcount++;
                    }
                    else if(nodeconnections[itemp][1] == p){
                        if(helpcount == q){
                            altid[q][0] = nodeindex[nodeconnections[itemp][0]][0];
                            altid[q][1] = nodeindex[nodeconnections[itemp][0]][1];
                            altid[q][2] = nodeindex[nodeconnections[itemp][0]][2];
                            altloc[q] = itemp;
                            break;
                        }
                        helpcount++;
                    }
                }
                //assign found point:)
                ti = altid[q][0];
                tj = altid[q][1];
                tk = altid[q][2];
                comppt = (sd->sB[ti][tj][tk]).nodepoint;
                //We now calculate the angle of the points
                if(*dim == 2){
                    //2D angle
                    double dx = comppt[0] - curpt[0];
                    double dy = comppt[1] - curpt[1];
                    angles[q] = atan2(dy,dx);
                    if(angles[q] < 0){
                        angles[q] = angles[q] + (2 * PI);
                    }
                }
                else{
                    //3D angle not implemented
                }
            }
            //Now we have each location & angle beween nodes
            //So we will sort the angles in a needed fashion
            //calculate the midpoint angles, and sort each node which falls between the angles
            sortAng(&angles,&altid,&altloc,&cid,t);
            double *newangles = malloc(cid * sizeof(double));
            for(int q = 0; q < cid; q++){
                //each new angle will fall between q & q + 1
                if(q != (cid - 1)){
                    newangles[q] = (angles[q] + angles[q + 1]) / 2;
                }
                else{
                    newangles[q] = ((angles[q] + angles[0]) / 2) + PI;
                    if(newangles[q] > 2*PI){
                        newangles[q] = newangles[q] - 2*PI;
                    }
                }
            }
            //We have calculated the midpoint angles; now we go though and assign nodes to splines
            double a1;
            double a2;
            for(int q = 0; q < cid; q++){
                //grab relavant angles
                if(q == 0){
                    a1 = newangles[cid - 1];
                    a2 = newangles[q];
                }
                else{
                    a1 = newangles[q - 1];
                    a2 = newangles[q];
                }
                //fprintf(stdout,"[%d,%d,%d] -> [%f,%f]\n",);
                for(int ni = 0; ni < *(sd->row); ni++){
                    for(int nj = 0; nj < *(sd->col); nj++){
                        for(int nk = 0; nk < *(sd->dep); nk++){
                            //check if node exists
                            if(*(sd->sB[ni][nj][nk]).leng != 0){
                                //Finally we will check if the node belongs to our current nodeid
                                if(((sd->sB[ni][nj][nk]).closeid1[0]==nodeindex[p][0])&&((sd->sB[ni][nj][nk]).closeid1[1]==nodeindex[p][1])&&((sd->sB[ni][nj][nk]).closeid1[2] == nodeindex[p][2])){
                                    //Only hitting relevant nodes now
                                    //so we can check angle according to id's
                                    if(ni==nodeindex[p][0]&&nj==nodeindex[p][1]&&nk==nodeindex[p][2]){
                                        //always add root nodes points :)
                                        for(int addpt = 0; addpt < *(sd->sB[ni][nj][nk]).leng; addpt++){
                                            findpoints[altloc[q]][comboindex[altloc[q]]] = (sd->sB[ni][nj][nk]).points[addpt];
                                            comboindex[altloc[q]]++;
                                        }
                                    }
                                    else{
                                        //Check angle if not rootnode, angle based upon index :)
                                        bool checkpass = false;
                                        double nangle = 0.;
                                        comppt = (sd->sB[ni][nj][nk]).points[0];
                                        if(*dim == 2){
                                            double dx = comppt[0] - curpt[0];
                                            double dy = comppt[1] - curpt[1];
                                            nangle = atan2(dy,dx);
                                            if(nangle < 0){
                                                nangle = nangle + 2 * PI;
                                            }
                                        }
                                        else{
                                        }
                                        //finally check for pass
                                        if(a1 < a2){
                                            if(nangle >= a1 && nangle <= a2){
                                                checkpass = true;
                                            }
                                        }
                                        else{
                                            if(nangle >= a1 || nangle <= a2){
                                                checkpass = true;
                                            }
                                        }
                                        //if passes add
                                        if(checkpass){
                                            for(int addpt = 0; addpt < *(sd->sB[ni][nj][nk]).leng; addpt++){
                                                findpoints[altloc[q]][comboindex[altloc[q]]] = (sd->sB[ni][nj][nk]).points[addpt];
                                                comboindex[altloc[q]]++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for(int q = 0; q < cid; q++){
                free(altid[q]);
            }
            free(altid);
            free(altloc);
            free(angles);
            free(newangles);
        }
    }
    free(idcount);
    free(sidcount);
    *pfindpoints = findpoints;
    *pnodeindex = nodeindex;
    *pnodeconnections = nodeconnections;
    *pcomboindex = comboindex;
    *sD = sd;
}
void pathNodes(struct skeleDensity **sD,int *dim,int **nodeid, int *nodeidcount,double ****pfindpoints,int **pcomboindex,int ***pnodeconnections,int ***pnodeindex,int *combocount,int *nicount,int *mxpt, double t){
    //First we Go through and Find which grid spaces have which id already tagged, and determine which nodes are touching which nodes,
    //This will work similar to the flood but stop when 
    struct skeleDensity *sd = *sD;
    double ***findpoints = *pfindpoints;
    int **nodeindex = *pnodeindex;
    int **nodeconnections = *pnodeconnections;
    int *comboindex = *pcomboindex;
    for(int i = 0; i < *nodeidcount; i++){
        int **visited = malloc((*(sd->pcount) + 1) * sizeof(int*));
        for(int j = 0; j < (*(sd->pcount) + 1); j++){
            visited[j] = calloc(3,sizeof(int));
        }
        visited[0][0] = nodeid[i][0];
        visited[0][1] = nodeid[i][1];
        visited[0][2] = nodeid[i][2];
        int visitcount = 0;
        floodDensity2(sD,dim,combocount,&nodeconnections,&nodeindex,nicount,&nodeid[i],&visited,&visitcount);
        for(int j = 0; j < *(sd->pcount); j++){
            free(visited[j]);
        }
        free(visited);
    }
    //With our second flood density, we now have the connections of the skeleton nodes ids
    //and thus we create paths, and assign nodes to connections
    findpoints = malloc(*combocount * sizeof(double**));
    comboindex = calloc(*combocount , sizeof(int));
    for(int i = 0; i < *combocount; i++){
        findpoints[i] = malloc(*mxpt * sizeof(double*));
    }
    sortConnections(sD,dim,combocount,&nodeconnections,&findpoints,&comboindex,&nodeindex,nicount,t);
    //output data
    //connection ID's
    char conname[80];
    sprintf (conname, "dat/connectionDat-%5.3f.dat", t);
    FILE * fpcon = fopen (conname, "w");
    for(int i = 0; i < *combocount;i++){
        fprintf(fpcon,"%d %d\n",nodeconnections[i][0],nodeconnections[i][1]);
    }
    fflush(fpcon);
    fclose(fpcon);
    //Connection Index:
    char indxname[80];
    sprintf (indxname, "dat/connectionidDat-%5.3f.dat", t);
    FILE * fpindx = fopen (indxname, "w");
    for(int i = 0; i < *nicount;i++){
        struct skeleBounds sb = (sd->sB[nodeindex[i][0]][nodeindex[i][1]][nodeindex[i][2]]);
        fprintf(fpindx,"%d %d %f %f\n",nodeindex[i][0],nodeindex[i][1],sb.nodepoint[0],sb.nodepoint[1]);
    }
    fflush(fpindx);
    fclose(fpindx);
    //reasign
    *pfindpoints = findpoints;
    *pnodeindex = nodeindex;
    *pnodeconnections = nodeconnections;
    *pcomboindex = comboindex;
    *sD = sd;
}
//calcs needed bezier properties for optimization at given t
double* calcBezierDC(int *n,double ***ppoints,double *t){
    double **points = *ppoints;
    if(*n == 1){
        double *retpoint = malloc(3*sizeof(double));
        retpoint[0] = points[0][0];
        retpoint[1] = points[0][1];
        retpoint[2] = points[0][2];
        return retpoint;
    }
    else{
        //malloc our new set of points we will calulate
        double **newpoints = malloc((*n - 1) * sizeof(double*));
        for(int i = 0; i < (*n - 1); i++){
            newpoints[i] = malloc(3 * sizeof(double));
            //De casteljau method
            newpoints[i][0] = (1 - *t) * points[i][0] + *t * points[i + 1][0];
            newpoints[i][1] = (1 - *t) * points[i][1] + *t * points[i + 1][1];
            newpoints[i][2] = (1 - *t) * points[i][2] + *t * points[i + 1][2];
        }
        int newn = *n - 1;
        //go one step deeper
        double *retpoint = calcBezierDC(&newn, &newpoints, t);
        for(int i = 0; i < (*n - 1); i++){
            free(newpoints[i]);
        }
        free(newpoints);
        return retpoint;
    }
}
double calcBezierErr(struct kdleaf *kdstruct,double **comppoints, int lengpoints, int *dim){
    struct kdleaf searchstruct = *kdstruct;
    double error = 0.;
    double *ignorepoint = malloc(3*sizeof(double));
    ignorepoint[0]=100000.;ignorepoint[1]=200000;ignorepoint[2]=400200;
    for(int i = 0; i < lengpoints; i++){
        double thiserror = 0.;
        int tleng = 1;
        double lowestdis = 0.;
        double *nearpoint = getNearest(comppoints[i],&searchstruct,&lengpoints,&ignorepoint,&tleng,&lowestdis);
        //position error 
        for(int j = 0; j < *dim + 1; j++){
            double dif = fabs(nearpoint[j] - comppoints[i][j]);
            thiserror += pow(dif,2);
        }
        error += thiserror;
    }
    error = error / lengpoints;
    free(ignorepoint);
    return error;
}
//fit spline itteratively
void findBestFit(double ***ppositioncoeff,double ***pradcoeff,double **findpoints,int comboindex,int *dim, int *n,double *tolerance,int splinediv,double inLam, double inDel){
    double **positioncoeff = *ppositioncoeff;
    double **radcoeff = *pradcoeff;
    //We have a few things we want to do here, firstly, we want to create a spline for the given section, which can be closest to our data points
    //so we create needed parameters
    double dt = 1/((double)splinediv);//NOTE t changes 0->1
    int maxitt = 100;
    double error = 1. + *tolerance;
    double error_last = 2. + *tolerance;
    double error_lastlast = 3. + *tolerance;
    double error_adjust = 2. + *tolerance;
    double lambda = inLam;
    double** deltas = malloc((*n - 1) * sizeof(double*));
    for(int i = 0; i < *n - 1; i++){
        deltas[i] = malloc((*dim + 1) * sizeof(double));
        for(int j = 0; j < *dim + 1; j++){
            //makes individual delta for each point and each 
            deltas[i][j] = inDel;
        }
    }
    //First we genrate our input values into a calculation spline
    double **calcspline = malloc((*n + 1) * sizeof(double*));
    for(int i = 0; i < (*n + 1); i++){
        calcspline[i] = malloc((*dim + 1) * sizeof(double));
        calcspline[i][0] = positioncoeff[i][0];
        calcspline[i][1] = positioncoeff[i][1];
        calcspline[i][2] = radcoeff[i][0];
        //fprintf(stdout,"init-spline:%d = [%f,%f,%f]\n",i,calcspline[i][0],calcspline[i][1],calcspline[i][2]);
    }
    int itt = 0;
    while(itt < maxitt && error > *tolerance){
        //Itterate through each middle point of the data at each dim
        int track = 0;
        //fprintf(stdout,"lambda:%f\n",lambda);
        for(int i = 1; i < (*n); i++){
            error_adjust = error;
            for(int j = 0; j < *dim + 1; j++){
                //define variables
                double temp = calcspline[i][j];
                double error_pos = 0.;
                double error_neg = 0.;
                //try positive adjustment
                calcspline[i][j] = temp + deltas[i - 1][j];
                double **ntpoints = malloc(splinediv * sizeof(double*));
                for(int p = 0; p < splinediv; p++){
                    double t = dt * (double)p;
                    int newn = *n + 1;
                    ntpoints[p] = malloc((*dim + 1)*sizeof(double));
                    double* newt = calcBezierDC(&newn,&calcspline,&t);
                    ntpoints[p][0] = newt[0];
                    ntpoints[p][1] = newt[1];
                    ntpoints[p][2] = newt[2];
                    free(newt);
                }
                struct kdleaf *kdstruct_petp = NULL;
                int splinecalc = splinediv - 1;
                CreateStructure(ntpoints,&kdstruct_petp,0,0,splinecalc,1);
                error_pos = calcBezierErr(kdstruct_petp,findpoints,comboindex,dim);
                kdDestroy(&kdstruct_petp);
                //Next try for negative adjustment
                calcspline[i][j] = temp - deltas[i - 1][j];
                for(int p = 0; p < splinediv; p++){
                    double t = dt * (double)p;
                    int newn = *n + 1;
                    double* newt = calcBezierDC(&newn,&calcspline,&t);
                    ntpoints[p][0] = newt[0];
                    ntpoints[p][1] = newt[1];
                    ntpoints[p][2] = newt[2];
                    free(newt);
                }
                struct kdleaf *kdstruct_petn = NULL;
                CreateStructure(ntpoints,&kdstruct_petn,0,0,splinecalc,1);
                error_neg = calcBezierErr(kdstruct_petn,findpoints,comboindex,dim);
                kdDestroy(&kdstruct_petn);
                //free points
                for(int k = 0; k < splinediv; k++){
                    free(ntpoints[k]);
                }
                free(ntpoints);
                //Determine next steps
                if(error_pos <= error_neg && error_pos < error_adjust){
                    //found better solution in the positive  
                    //we mark a point for decreasing scale later too
                    error_adjust = error_pos;
                    calcspline[i][j] = temp + deltas[i - 1][j] * lambda;
                    track++;
                    //because sucess we decrease delta for this parameter as we get closer to the solution
                    deltas[i - 1][j] = deltas[i - 1][j] * 0.1;
                    //fprintf(stdout,"(%d-%d) => new pos [%f,%f,%f]\n",i-1,j,calcspline[i][0],calcspline[i][1],calcspline[i][2]);
                }
                else if(error_neg < error_pos && error_neg < error_adjust){
                    //found better solution in the negative 
                    //we mark a point for decreasing scale later too
                    error_adjust = error_neg;
                    calcspline[i][j] = temp - deltas[i - 1][j] * lambda;
                    track++;
                    //because sucess we decrease delta for this parameter as we get closer to the solution
                    deltas[i - 1][j] = deltas[i - 1][j] * 0.1;
                    //fprintf(stdout,"(%d-%d) => new  pos [%f,%f,%f]\n",i-1,j,calcspline[i][0],calcspline[i][1],calcspline[i][2]);
                }
                else{
                    //neither was a better shot, so now we reset  
                    //we also keep delta the same
                    calcspline[i][j] = temp;
                    //fprintf(stdout,"(%d-%d) => miss pos\n",i-1,j);
                }
            }
        }
        if(track > 0){
            //finally at the end if our spline is closer we decrease lambda
            lambda = lambda * 0.1;
            error = error_adjust;
            if(error_last == error && error_lastlast == error){
                //breakout as we are stuck
                break;
            }
        }
        else{
            //we never did better than original error, so we increase lambda
            if(error_last == error && error_lastlast == error){
                //breakout as we are stuck
                break;
            }
            else{
                lambda = lambda * 10;
            }
        }
        printf("itt:%d|l:%f|e:%f|el:%f|ell:%f \n",itt,lambda,error,error_last,error_lastlast);
        error_lastlast = error_last;
        error_last = error;
        itt++;
    }
    printf("finished at:%d\n",itt);
    for(int i = 0; i < (*n + 1); i++){
        positioncoeff[i][0] = calcspline[i][0];
        positioncoeff[i][1] = calcspline[i][1];
        radcoeff[i][0] = calcspline[i][2];
        //fprintf(stdout,"post-spline:%d = [%f,%f,%f]\n",i,calcspline[i][0],calcspline[i][1],calcspline[i][2]);
    }
    //cleanup
    for(int i = 0; i < (*n + 1); i++){
        free(calcspline[i]);
    }
    free(calcspline);
    for(int i = 0; i < *n - 1; i++){
        free(deltas[i]);
    }
    free(deltas);
    //force reassign for unbreakable update
    *ppositioncoeff = positioncoeff; 
    *pradcoeff = radcoeff;      
}
double getDetN(double ***pA,int *n){
    double **A = *pA;
    double ret = 0.;
    double factor = 1.;
    if(*n == 1){
        return **A;
    }
    for(int i = 0; i < *n; i++){
        //Malloc our new A
        double **newA = malloc((*n - 1) * sizeof(double*));
        for(int j = 0; j < (*n - 1); j++){
            newA[j] = malloc((*n - 1) * sizeof(double));
        }
        for(int j = 1; j < *n;j++){
            for(int k = 0; k < *n; k++){
                if(k==i) continue;
                newA[j-1][k<i?k:(k-1)]=A[j][k];
            }
        }
        int newn = *n - 1;
        ret += factor*A[0][i]*getDetN(&newA,&newn);
        factor *= -1.0;
        for(int j = 0; j < *n - 1; j++){
            free(newA[j]);
        }
        free(newA);
    }
    return ret;
}
//factorial
double getFactorial(int num){
    int ret = 1;
    for(int i = 1; i <= num; i++){
        ret *= i;
    }
    return (double)ret;
}
////carrier data :)
//struct bezData{
//    double **points;
//    double **Bval;
//    int *npoints;
//    int *cpoints;
//    int *dims;
//};
////alloc and dealloc functions for bezData
//void createbezData(struct bezData *bd,double **datapoints,int *datalength,double **bvalues,int *dim,int *degree){
//    //creating data struct
//    struct bezData newbez = malloc(sizeof(struct bezData));
//    newbez.points = datapoints;
//    newbez.Bval = bvalues;
//    newbez.dims = dim;
//    newbez.npoints = datalength;
//    newbez.cpoints = degree;
//    bd = &newbez;
//}
//void destroybezData(struct bezData *bd){
//    //destroying data struct
//    free(bd);
//}
//model functions :)
//double B_obj_func(const gsl_vector *x, void *params) {
//    struct bezData *pbd = (struct bezData*)params;
//    double **B = (*pbd).Bval;
//    double **data = (*pbd).points;
//    double dev = 0.;
//    int divmat = *(pbd->npoints) / *(pbd->dims);
//    for(int i = 0; i < *(pbd->npoints); i++){
//        //itterate through each point :0
//        double dd = 0.;
//        for(int j = 0; j < *(pbd->dims);j++){
//            double error = 0.;
//            double calcbez = 0.;
//            for(int k = 0; k < *(pbd->cpoints); k++){
//                //calc the sum bez point for comparison on dim
//                calcbez += B[i][k] * gsl_vector_get(x,k+j*divmat);
//            }
//            error = calcbez - data[i][j];
//            dd += pow(error,2);
//        }
//        dev += dd;
//    }
//    //total error = sqrt (1/N * sum errors square)
//    dev = sqrt( ( 1. / (double)(*(pbd->cpoints))) * (dev));
//    return dev;
//}
//void getCoeffBNLC(double **datapoints,int lengdata,double **B,int *dim,int n){
//    //solve bezier for non-linear constrained problem
//    //using multimin optimzation function
//    //Assing struct
//    struct bezData bd;
//    createbezData(&bd,datapoints,&lengdata,B,dim,&n);
//
//    //setup solver
//    gsl_multimin_function min_func;
//    min_func.n = *dim * n;
//    min_func.f = B_obj_func;
//    min_func.params = &bd;
//     
//    //free/destroy
//    destroybezData(&bd);
//}
void getCoeffBL(int row, int col, double ***pa, double **px,double **pp,int dim){
    //solve bezier for linear non constrained problem
    //finally we solve the equation A*P=X
    //alloc for gsl
    double **a = *pa;
    double *x = *px;
    double *p = *pp; 
    gsl_matrix *A = gsl_matrix_alloc(row,col); 
    gsl_vector *X = gsl_vector_alloc(row); 
    gsl_vector *P = gsl_vector_alloc(col);
    gsl_matrix *cov = gsl_matrix_alloc(row,row);
    int divmat = col / dim;
    printf("div=%d\n",divmat);
    for(int i = 0; i < dim; i++){
        gsl_vector_set(P,i*divmat,x[i*divmat]);
        gsl_vector_set(P,(divmat-1)+i*divmat,x[(divmat-1)+i*divmat]);
    }
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            gsl_matrix_set(A,i,j,a[i][j]);
        }
        gsl_vector_set(X,i,x[i]);
    }
    double chisq = 0.;
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(row, col);
    gsl_multifit_linear(A, X, P, cov, &chisq, work);
    gsl_multifit_linear_free(work);
    for(int i = 0; i < col; i++){
        (p[i]) = gsl_vector_get(P,i);
    }
    //free up values
    gsl_matrix_free(A);
    gsl_vector_free(X);
    gsl_vector_free(P);
    gsl_matrix_free(cov);
    *pa = a;
    *px = x;
    *pp = p;
}
void getCoeffBLW(int row, int col, double ***pa, double **px,double **pp,int dim){
    //solve bezier for linear non constrained *With Weights* problem
    //finally we solve the equation A*P=X
    //alloc for gsl
    double **a = *pa;
    double *x = *px;
    double *p = *pp; 
    gsl_matrix *A = gsl_matrix_alloc(row,col); 
    gsl_vector *X = gsl_vector_alloc(row); 
    gsl_vector *P = gsl_vector_alloc(col);
    gsl_vector *W = gsl_vector_alloc(row);
    gsl_matrix *cov = gsl_matrix_alloc(row,row);
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            gsl_matrix_set(A,i,j,a[i][j]);
        }
        gsl_vector_set(X,i,x[i]);
        if(i % (row / dim) == 0){
            gsl_vector_set(W,i,2);
        } 
        else if(i % (row / dim) == (row/dim) - 1){
            gsl_vector_set(W,i,2);
        }
        else{
            gsl_vector_set(W,i,0.1);
        }
    }
    double chisq = 0.;
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(row, col);
    gsl_multifit_wlinear(A, W, X, P, cov, &chisq, work);
    gsl_multifit_linear_free(work);
    for(int i = 0; i < col; i++){
        (p[i]) = gsl_vector_get(P,i);
    }
    //free up values
    gsl_matrix_free(A);
    gsl_vector_free(X);
    gsl_vector_free(P);
    gsl_vector_free(W);
    gsl_matrix_free(cov);
    *pa = a;
    *px = x;
    *pp = p;
}
void getCoeffBLWBC(int row, int col, double ***pa, double **px,double **pp,int dim,double **applyGradient,int aGradlen){
    //solve bezier for linear non constrained *With Weights* problem
    //finally we solve the equation A*P=X
    //alloc for gsl
    double **a = *pa;
    double *x = *px;
    double *p = *pp; 
    gsl_matrix *A = gsl_matrix_alloc(row,col); 
    gsl_vector *X = gsl_vector_alloc(row); 
    gsl_vector *P = gsl_vector_alloc(col);
    gsl_vector *W = gsl_vector_alloc(row);
    gsl_matrix *cov = gsl_matrix_alloc(row,row);
    int rd = row/dim;
    int cd = col/dim;
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            gsl_matrix_set(A,i,j,a[i][j]);
        }
        gsl_vector_set(X,i,x[i]);
        if(i % rd == 0){//at 0 positions
            gsl_vector_set(W,i,1);
        } 
        else if(i % rd == rd - 1){//at n positions
            gsl_vector_set(W,i,1);
        }
        else{
            gsl_vector_set(W,i,0.1);
        }
    }
    //Next we sort and apply boundary conditions
    if(aGradlen == 1){
        if(applyGradient[0][0] != 0. && applyGradient[1][0] == 0.){
            //if resctriction to start of spline
            for(int i = 0; i < dim; i++){
                //goes through each dimension
                for(int j = 0; j < rd; j++){
                    for(int k = 0; k < cd; k++){
                        if(k == 0){
                            gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)+gsl_matrix_get(A,j+i*rd,(k + 1)+i*cd));//sets B0 to be B0+B1
                        }
                        else if(k == 1){
                            gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)*applyGradient[0][i]);//sets B1 to to be B1*d
                        }
                        //Our system is now in the form of (B0*P0+B1*(P1=P0+TP1*der)+...+BnPn) thus all P0's now haeve a B0&B1 componenet
                        //P1 converts to tempPoint TP1 and der represents the direction vector defined for this point
                        //After fitting we must convert back into usable points using this P1=P0+TP1*der; 
                        //TP1 essentially will be the fitting distance along the gradient :)
                    }
                }
            }
            printf("sucess on gradient convert\n");
        }
        else if(applyGradient[1][0] != 0. && applyGradient[0][0] == 0.){
            //resctriciton on end of spline
            for(int i = 0; i < dim; i++){
                //goes through each dimension
                for(int j = 0; j < rd; j++){
                    for(int k = 0; k < cd; k++){
                        if(k == cd - 1){
                            gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)+gsl_matrix_get(A,j+i*rd,(k - 1)+i*cd));//sets Bn to be Bn+Bn-1
                        }
                        else if(k == cd - 2){
                            gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)*applyGradient[1][i]);//sets Bn-1 to be Bn-1*d
                        }
                    }
                }
            }
            printf("sucess on gradient convert 2 [error]\n");
        }
        else{
            printf("error no applicable gradient error\n");
        }
    }
    else{
        //we have a gradient on both vectors ;)
        //resctriciton on end of spline
        for(int i = 0; i < dim; i++){
            //goes through each dimension
            for(int j = 0; j < rd; j++){
                for(int k = 0; k < cd; k++){
                    if(k == 0){
                        gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)+gsl_matrix_get(A,j+i*rd,(k + 1)+i*cd));//sets B0 to be B0+B1
                    }
                    else if(k == 1){
                        gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)*applyGradient[0][i]);//sets B1 to to be B1*d
                    }
                    else if(k == cd - 1){
                        gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)+gsl_matrix_get(A,j+i*rd,(k - 1)+i*cd));//sets Bn to be Bn+Bn-1
                    }
                    else if(k == cd - 2){
                        gsl_matrix_set(A,j+i*rd,k+i*cd,gsl_matrix_get(A,j+i*rd,k+i*cd)*applyGradient[1][i]);//sets Bn-1 to be Bn-1*d
                    }
                }
            }
        }
        printf("sucess on double gradient convert\n");
    }
    double chisq = 0.;
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(row, col);
    gsl_multifit_wlinear(A, W, X, P, cov, &chisq, work);
    gsl_multifit_linear_free(work);
    //finally convert back to our control points
    if(aGradlen == 1){
        if(applyGradient[0][0] != 0. && applyGradient[1][0] == 0.){
            //if resctriction to start of spline
            for(int i = 0; i < dim; i++){
                //goes through each dimension
                for(int k = 0; k < cd; k++){
                    if(k == 1){
                        gsl_vector_set(P,k+i*cd,gsl_vector_get(P,(k-1)+i*cd) + gsl_vector_get(P,k+i*cd) * applyGradient[0][i]);//sets P1 to P1 = P0 + TP1 * d 
                    }
                }
            }
            printf("sucess on gradient convert  back\n");
        }
        else if(applyGradient[1][0] != 0. && applyGradient[0][0] == 0.){
            //resctriciton on end of spline
            for(int i = 0; i < dim; i++){
                //goes through each dimension
                for(int k = 0; k < cd; k++){
                    if(k == cd - 2){
                        gsl_vector_set(P,k+i*cd,gsl_vector_get(P,(k+1)+i*cd) + gsl_vector_get(P,k+i*cd) * applyGradient[1][i]);//sets Pn-1 to Pn-1 = Pn + TPn-1 * d 
                    }
                }
            }
            printf("sucess on gradient convert back\n");
        }
        else{
            printf("error no applicable gradient error\n");
        }
    }
    else{
        //we have a gradient on both vectors ;)
        //resctriciton on end of spline
        for(int i = 0; i < dim; i++){
            //goes through each dimension
            for(int k = 0; k < cd; k++){
                if(k == 1){
                    gsl_vector_set(P,k+i*cd,gsl_vector_get(P,(k-1)+i*cd) + gsl_vector_get(P,k+i*cd) * applyGradient[0][i]);//sets P1 to P1 = P0 + TP1 * d 
                }
                else if(k == cd - 2){
                    gsl_vector_set(P,k+i*cd,gsl_vector_get(P,(k+1)+i*cd) + gsl_vector_get(P,k+i*cd) * applyGradient[1][i]);//sets Pn-1 to Pn-1 = Pn + TPn-1 * d 
                }
            }
        }
        printf("sucess on double gradient convert back\n");
    }

    for(int i = 0; i < col; i++){
        (p[i]) = gsl_vector_get(P,i);
    }
    //free up values
    gsl_matrix_free(A);
    gsl_vector_free(X);
    gsl_vector_free(P);
    gsl_vector_free(W);
    gsl_matrix_free(cov);
    *pa = a;
    *px = x;
    *pp = p;
}
//Least square fitting splines
void findBestFit2(double ***ppositioncoeff,double ***pradcoeff,double **findpoints,int comboindex,int *dim, int *n,double *tolerance){
    double **positioncoeff = *ppositioncoeff;
    double **radcoeff = *pradcoeff;
    //we have B matrix which collects all of our desired B values for each datapoint
    double **B = malloc((comboindex)*sizeof(double*));
    for(int i = 0; i < comboindex; i++){
        B[i] = calloc((*n + 1),sizeof(double));
    }
    //construct t for finding values
    double *t = malloc(comboindex * sizeof(double));
    int calcskip = 0;
    for(int i = 0; i < comboindex; i++){
        //for each data point try and guess the t value of the data based on linear sampling
        //allows for least square method to match up best guess spline points to data better based on flucuations
        //int denomt = *dim + 1;
        //double tx = (positioncoeff[0][0] - findpoints[i][0]) / (positioncoeff[0][0] - positioncoeff[*n][0]);
        //if(tx > 1.){
        //    tx = 0.99;
        //}
        //else if(tx < 0.){
        //    tx = 0.01;
        //}
        //double ty = (positioncoeff[0][1] - findpoints[i][1]) / (positioncoeff[0][1] - positioncoeff[*n][1]);
        //if(ty > 1.){
        //    ty = 0.99;
        //}
        //else if(ty < 0.){
        //    ty = 0.01;
        //}
        //double tr = (radcoeff[0][0] - findpoints[i][2]) / (radcoeff[0][0] - radcoeff[*n][0]);
        //if(tr > 1.){
        //    tr = 0.99;
        //}
        //else if(tr < 0.){
        //    tr = 0.01;
        //}
        if(*dim == 2){
            t[i] = (double)(i) / (double)comboindex;
            //fprintf(stdout,"got t-%d/%d: %f\n",i,comboindex,t[i]);
            //if the calced t value is a bad take or past the endpoint average, we will 
            //if(denomt == 0){
            //    //all values have no linear t value :(
            //    t[i] = -1.;
            //    calcskip++;
            //}
            //else{
            //    t[i] = (tx + ty + tr) / denomt;
      
            //if(t[i] < 0.01 || t[i] > 0.99){
            //    //    t[i] = -1.;
            //    //    calcskip++;
            //    //}
            //}
            //fprintf(stdout,"gott t-%d: %f = (%f + %f + %f) / %d\n",i,t[i],(positioncoeff[0][0] - findpoints[i][0]) / (positioncoeff[0][0] - positioncoeff[*n][0]),(positioncoeff[0][1] - findpoints[i][1]) / (positioncoeff[0][1] - positioncoeff[*n][1]),(radcoeff[0][0] - findpoints[i][2]) / (radcoeff[0][0] - radcoeff[*n][0]),denomt);
        }
        else{
            printf("error not implemented-3D :(((( \n");
        }
    }
    t[0] = 0.;//apply boundary equation to system
    t[comboindex-1] = 1.;//initzalize correct start & end chord length
    //Next we calc our B function:)
    for(int i = 0; i < comboindex; i++){
        //blocks off bad t's from being calc
        //if(t[i] != -1.){
        double tt = 1. - t[i];
        //Find B values at current points approx
        for(int j = 0; j < (*n + 1); j++){
            B[i][j] = pow(tt,(*n) - j) * pow(t[i],(j)) * (getFactorial(*n) / (getFactorial(j) * getFactorial((*n) - j)));
            //fprintf(stdout,"calcB%d(%f) => %f = %f * %f * %f\n",j,t[i],B[i][j],pow(tt,(*n) - j),pow(t[i],(j)),(getFactorial(*n) / (getFactorial(j) * getFactorial((*n) - j))));
        }
        //}
    }
    //solve system
    int newci = comboindex - calcskip;
    //Solve all 3 Dim + extra at once for no BC
    int solscale = *dim + 1;
    double **calcB = malloc((solscale * newci)*sizeof(double*));
    double *calcD = calloc(solscale * newci,sizeof(double));
    double *calcP = calloc(solscale * (*n + 1),sizeof(double));
    for(int i = 0; i < solscale*newci; i++){
        calcB[i] = calloc(solscale * (*n + 1),sizeof(double));
    }
    for(int q = 0; q < solscale; q++){
        for(int i = 0; i < newci; i++){
            for(int j = 0; j < *n + 1; j++){
                calcB[i+q*newci][j+q*(*n + 1)] = B[i][j];
            }
            if(i == 0){
                if(q < *dim){
                    calcD[i+q*newci] = positioncoeff[0][q];
                }
                else{
                    calcD[i+q*newci] = radcoeff[0][q-*dim];
                }
            }
            else if(i == newci - 1){
                if(q < *dim){
                    calcD[i+q*newci] = positioncoeff[*n][q];
                }
                else{
                    calcD[i+q*newci] = radcoeff[*n][q-*dim];
                }
            }
            else{
              calcD[i+q*newci] = findpoints[i][q];
            }
        }
    }
    //before pushing to the solution we consider adding in boundary conditions
    getCoeffBLW(solscale*newci,solscale*(*n+1),&calcB,&calcD,&calcP,solscale);
    //if(aGradlen > 0){
    //    //note we adjust solution variables inside of the function, we get output useable variables in the same format
    //    getCoeffBLWBC(solscale*newci,solscale*(*n+1),&calcB,&calcD,&calcP,solscale,applyGradient,aGradlen);
    //}
    //else{
    //}
    for(int i = 0; i < solscale; i++){
        for(int j = 1; j < *n; j++){
            //fprintf(stdout,"P%d:%d => %f\n",j,i,calcP[j + i * (*n + 1)]);
            if(i < *dim){
                positioncoeff[j][i] = calcP[j + i * (*n + 1)]; 
            }
            else{
                radcoeff[j][i-*dim] = calcP[j + i * (*n + 1)]; 
            }
        } 
    }
    for(int j = 0; j < solscale*newci; j++){
        free(calcB[j]);
    }
    free(calcB);
    free(calcD);
    free(calcP);
    //for(int i = 0; i < *dim + 1; i++){
    //    //Next we create & solve a linear system along each dimension
    //    double **calcB = malloc((newci) * sizeof(double*));
    //    double *calcD = calloc((newci) , sizeof(double));
    //    double *calcP = calloc((*n + 1) , sizeof(double));
    //    //fprintf(stdout,"\n");
    //    for(int j = 0; j < *n + 1; j++){
    //        //fprintf(stdout,"calcP%d=%f\n",j,calcP[j]);
    //    }
    //    //fprintf(stdout,"\n");
    //    int j = 0;
    //    for(int ittj = 0; ittj < comboindex; ittj++){
    //        //blocks off bad t's from being calc
    //        if(t[ittj] != -1.){
    //            calcB[j] = calloc(*n + 1,sizeof(double));
    //            for(int q = 0; q < (*n + 1); q++){
    //                //assign calcB values
    //                calcB[j][q] = B[ittj][q];
    //            }
    //            calcD[j] = findpoints[ittj][i];
    //            j++;
    //        }
    //    }
    //    //finally set BC's on calcP :)
    //    if(i < *dim){
    //        calcP[0]  = positioncoeff[0][i];
    //        calcP[*n] = positioncoeff[*n][i];
    //    }
    //    else{
    //        calcP[0]  = radcoeff[0][0];
    //        calcP[*n] = radcoeff[*n][0];
    //    }
    //    //int midpoint = (int)ceil(newci/2);
    //    //fprintf(stdout,"Starting Calc on Dim %d\n",i);
    //    //for(int j = 0; j < newci; j++){
    //    //    //print B portion of matrix
    //    //    fprintf(stdout,"[ ");
    //    //    for(int k = 0; k < (*n + 1); k++){
    //    //        if(k != *n){
    //    //            fprintf(stdout,"%f , ",calcB[j][k]);
    //    //        }
    //    //        else{
    //    //            fprintf(stdout,"%f ",calcB[j][k]);
    //    //        }
    //    //    }
    //    //    fprintf(stdout,"] ");
    //    //    //print P's and spaces
    //    //    if(j == midpoint && midpoint < (*n + 1)){
    //    //        fprintf(stdout," [ P%d:%d ]  =====>  ",midpoint,i);
    //    //    }
    //    //    else if( j == midpoint){
    //    //        fprintf(stdout,"          =====>     ");
    //    //    }
    //    //    else if(j < (*n + 1)){
    //    //        if(j == 0){
    //    //            fprintf(stdout," [ P%d:%d %f ]            ",j,i,calcP[0]);
    //    //        }
    //    //        else if(j == *n){
    //    //            fprintf(stdout," [ P%d:%d %f ]            ",j,i,calcP[*n]);
    //    //        }
    //    //        else{
    //    //            fprintf(stdout," [ P%d:%d ]            ",j,i);
    //    //        }
    //    //    }
    //    //    else{
    //    //        fprintf(stdout,"                     ");
    //    //    }
    //    //    //print D matrix
    //    //    fprintf(stdout,"[ %f ] \n",calcD[j]);
    //    //}
    //    //now we have defined both calcD and calcP so we will use gausian elimation to find the solutions
    //    if(newci> *n+1){
    //        getCoeffB(newci,*n+1,&calcB,&calcD,&calcP);
    //        for(int j = 1; j < *n; j++){
    //            fprintf(stdout,"P%d:%d => %f\n",j,i,calcP[j]);
    //            if(i < *dim){
    //                positioncoeff[j][i] = calcP[j]; 
    //            }
    //            else{
    //                radcoeff[j][i-*dim] = calcP[j]; 
    //            }
    //        }
    //    }
    //    for(int j = 0; j < newci; j++){
    //        free(calcB[j]);
    //    }
    //    free(calcB);
    //    free(calcD);
    //    free(calcP);
    //}
    //for(int i = 0; i < *n - 1; i++){
    //    //finally we will go though each point we want to approximate and assign it
    //    //i represents each output point
    //    double num;
    //    for(int j = 0; j < *dim; j++){
    //        num = 0.;
    //        for(int k = 0; k < *n - 1; k++){
    //            num += (A[k][i]) * coldims[k][j];
    //        }
    //        //itterate through each dim
    //        positioncoeff[i+1][j] = num / denom;
    //    }
    //    num = 0.;
    //    for(int k = 0; k < *n - 1; k++){
    //        num += (A[k][i]) * coldims[k][*dim];
    //    }
    //    radcoeff[i+1][0] = num / denom;
    //}
    //Free
    free(t);
    free(B);
    *ppositioncoeff = positioncoeff; 
    *pradcoeff = radcoeff;      
}
void getSplineCoeff(double ***ppositioncoeff,double ***pradcoeff,double **pnode0,double **pnode1,int *dim,int *n,double ***pfindpoints,int *pcomboindex,double *tolerance){
    double **positioncoeff = *ppositioncoeff;
    double **radcoeff = *pradcoeff;
    double *node0 = *pnode0;
    double *node1 = *pnode1;
    double **findpoints = *pfindpoints;
    int comboindex = *pcomboindex;
    //To create our bezier curve, we will first generate needed controlpoints
    if(*n == 1){
        //if n == 1, then we dont need to generate controlpoints
        positioncoeff[0][0] = node0[0];//since linear starting node
        positioncoeff[0][1] = node0[1];
        radcoeff[0][0] = node0[2];
        positioncoeff[1][0] = node1[0];//ending node
        positioncoeff[1][1] = node1[1];
        radcoeff[1][0] = node1[2];
    }
    else{
        //if n is bigger than we need to find new bestfit
        //assign start & end
        positioncoeff[0][0] = node0[0];//starting node
        positioncoeff[0][1] = node0[1];
        radcoeff[0][0] = node0[2];
        positioncoeff[*n][0] = node1[0];//ending node
        positioncoeff[*n][1] = node1[1];
        radcoeff[*n][0] = node1[2];
        //next we create initial points for our spline, we inizalize each section as 'linear'
        double dx = (node1[0] - node0[0])/((double)*n);
        double dy = (node1[1] - node0[1])/((double)*n);
        double dr = (node1[2] - node0[2])/((double)*n);
        for(int i = 1; i < *n; i++){
            positioncoeff[i][0] = positioncoeff[i - 1][0] + dx;
            positioncoeff[i][1] = positioncoeff[i - 1][1] + dy;
            radcoeff[i][0] = radcoeff[i - 1][0] + dr;
        }
        //Only Output spline if enough points, otherwise we change
        if(comboindex > *n){
            //findBestFit(&positioncoeff,&radcoeff,findpoints,comboindex,dim,n,tolerance,30,inLam,inDel);
            printf("Spline before: ");
            for(int i = 0; i < *n + 1; i++){
                printf("[%f,%f,%f]",positioncoeff[i][0],positioncoeff[i][1],radcoeff[i][0]);
            }
            printf("\n");
            findBestFit2(&positioncoeff,&radcoeff,findpoints,comboindex,dim,n,tolerance);
            printf("Spline after: ");
            for(int i = 0; i < *n + 1; i++){
                printf("[%f,%f,%f]",positioncoeff[i][0],positioncoeff[i][1],radcoeff[i][0]);
            }
            printf("\n\n");
        }
        else{
            printf("error too little input points defaulting spline to linear/input\n");
            //if(comboindex == 1){
            //    for(int i = 1; i < *n; i++){
            //        positioncoeff[i][0] = findpoints[0][0];
            //        positioncoeff[i][1] = findpoints[0][1];
            //        radcoeff[i][0] = findpoints[0][2];
            //    }
            //}
            //else{
            //    int mid = (int)floor((*n + 1) / 2);
            //    for(int i = 1; i < *n; i++){
            //        if(i < mid){
            //            positioncoeff[i][0] = findpoints[0][0];
            //            positioncoeff[i][1] = findpoints[0][1];
            //            radcoeff[i][0] = findpoints[0][2];
            //        }
            //        else{
            //            positioncoeff[i][0] = findpoints[comboindex - 1][0];
            //            positioncoeff[i][1] = findpoints[comboindex - 1][1];
            //            radcoeff[i][0] = findpoints[comboindex - 1][2];
            //        }
            //    }

            //}
        }
    }
    *ppositioncoeff = positioncoeff; 
    *pradcoeff = radcoeff;      
}
//Splining method, for merging unessicary nodes, and creating a spline
void makeSpline(struct skeleDensity **sD,int *dim,double *tolerance,int *length,double t,int *n){
    //Here the only node points are Skeleton points ends, we move from these
    //along the region of interest 
    struct skeleDensity *sDmain = *sD;
    //for holding our locations & amounts of endpoints
    int ncount = *sDmain->ncount;
    int **nodeid = malloc(ncount * sizeof(int*));
    for(int i = 0; i < ncount; i++){
        nodeid[i] = calloc(3,sizeof(int));
    }
    //first we collect the locations of each node
    int localcount = 0;
    int stackcount = 0;//allocates max needed space in our stack
    //make adjustments if needed
    printf("nodepoints b4 = %d\n",ncount);
    if(*(sDmain->row) > 5 || *(sDmain->col) > 5 || *(sDmain->dep) > 5){
        reduceLocalMax(sD,dim);
    }
    printf("nodepoints af = %d\n",*sDmain->ncount);
    for(int i = 0; i < *sDmain->row; i++){
        for(int j = 0; j < *sDmain->col; j++){
            for(int k = 0; k < *sDmain->dep; k++){
                if(*(sDmain->sB[i][j][k]).hasNode){// && *(sDmain->sB[i][j][k]).smode  == 2){
                    nodeid[localcount][0] = i; 
                    nodeid[localcount][1] = j; 
                    nodeid[localcount][2] = k;
                    localcount++;
                }
                if(*(sDmain->sB[i][j][k]).leng > 0){
                    stackcount++;
                }
            }
        }
    }
    //Here we begin our distance calculations from each endpoint
    int **visitstack = malloc(stackcount * stackcount * sizeof(int*));
    for(int i = 0; i < stackcount * stackcount; i++){
        visitstack[i] = calloc(3,sizeof(int));
    }
    int *visitcount = calloc(1,sizeof(int));
    int idcount = localcount-1;
    int **ids = malloc(idcount * sizeof(int*));
    for(int i = 0; i < idcount; i++){
        ids[i] = calloc(3,sizeof(int));
    }
    for(int i = 0; i < localcount; i++){
        *visitcount = 0;
        int countnow = 0;
        for(int j = 0; j < localcount; j++){
            if(j != i){
                ids[countnow][0] = nodeid[j][0];
                ids[countnow][1] = nodeid[j][1];
                ids[countnow][2] = nodeid[j][2];
                countnow++;
            }
        }
        floodDensity(sD,dim,nodeid[i],0,&visitstack,visitcount,ids,&idcount,tolerance);
        if(i != idcount){
            resetFloodStack(visitstack,stackcount * stackcount,ids,idcount);
        }
    }
    //We have flooded our tree with endpoint calculations,
    //Next we want to find the paths between the points :)
    double ***findpoints = NULL;//A collection of points which will be combined with the node id's[[pts1],[pts2],[pts3]]
    int **nodeconnections = NULL;//A collection of nodes tied to points [[i1,i2],[i1,i2],[i3,i4],..]
    int **nodeindex = NULL;//A index for the nodes [[32,54,0],[2,6,8],..]
    int *comboindex = NULL;//current index of each individual combo
    int combocount = 0;
    int nicount = 0;
    pathNodes(sD,dim,nodeid,&localcount,&findpoints,&comboindex,&nodeconnections,&nodeindex,&combocount,&nicount,length,t);
    //Now weve assigned nodes & collected splines approx points, we will combine branches = 2 into larger sections
    printf("combo count => %d\n",combocount);
    if(combocount > 1){
        thinNodePoint(sD,dim,&findpoints,&comboindex,&nodeconnections,&nodeindex,&combocount,&nicount,length,t);
    }
    //Next we will go through and calculate our approximations
    //first we define some things we want
    //given we want order 'n' we allocate accordingly
    
    //allocation for coeffs
    double ***positioncoeff = malloc(combocount*sizeof(double**));
    double ***radcoeff = malloc(combocount*sizeof(double**));//NOTE rad coeff is set to be build like position coeff to allow approximation of other var
    double **node0 = malloc(combocount*sizeof(double*));
    double **node1 = malloc(combocount*sizeof(double*));
    double mindis = *sDmain->xmax - *sDmain->xmin;
    double mindisy = *sDmain->ymax - *sDmain->ymin;
    if(mindisy < mindis){
        mindis = mindisy;
    }
    if(*dim == 3){
        double mindisz = *sDmain->ymax - *sDmain->ymin;
        if(mindisz < mindis){
            mindis = mindisz;
        }
    }
    int *newn = calloc(combocount,sizeof(int));
    //double **markGradient = malloc(localcount * sizeof(double*));//markgradient holds the index positio's node gradient :) 
    //int *hasGradient = calloc(localcount , sizeof(int));//holds 0 for not calculated; 1 for 2 connections or wanted gradient; -1 for 1 or more than 2 connections
    //for(int q = 0; q  < localcount; q++){
    //    fprintf(stdout,"markGrad alloc %d-%d\n",q,*dim + 1);
    //    markGradient[q] = calloc((*dim + 1) , sizeof(double));
    //}
    for(int q = 0; q < combocount; q++){
        //add already calculated gradients into solution steps :)
        //note keeps one and then two index, fills other with zero; 
        //double **applyGradient = malloc(2*sizeof(double*));
        //int aGradlen = 0;
        //if(hasGradient[nodeconnections[q][0]] == 1){
        //    applyGradient[0] = markGradient[nodeconnections[q][0]];
        //    aGradlen++;
        //}
        //else{
        //    applyGradient[0] = markGradient[nodeconnections[q][0]];
        //}
        //if(hasGradient[nodeconnections[q][1]] == 1){
        //    applyGradient[1] = markGradient[nodeconnections[q][1]];
        //    aGradlen++;
        //}
        //else{
        //    applyGradient[1] = markGradient[nodeconnections[q][1]];
        //}

        //Now we have the space for our current coefficients so next we grab the relevant points and calculate them
        int extra = 1;
        int i0=nodeindex[nodeconnections[q][0]][0],j0=nodeindex[nodeconnections[q][0]][1],k0=nodeindex[nodeconnections[q][0]][2];
        int i1=nodeindex[nodeconnections[q][1]][0],j1=nodeindex[nodeconnections[q][1]][1],k1=nodeindex[nodeconnections[q][1]][2];
        node0[q] = malloc((*dim + extra)*sizeof(double));
        node1[q] = malloc((*dim + extra)*sizeof(double));
        for(int p = 0; p < *dim + extra; p++){
            node0[q][p] = (sDmain->sB[i0][j0][k0]).nodepoint[p];
            node1[q][p] = (sDmain->sB[i1][j1][k1]).nodepoint[p];
        }
        //Next we rank n based on knowndistance, our smallest dimensional parameter
        double dist = getDistance(node0[q],node1[q]);//distance between known points
        //fprintf(stdout,"mindis=%f; dis=%f\n",mindis,dist);
        //next we find min distance of sD
        newn[q] = (int)ceil(*n * dist / mindis);
        if(newn[q] > *n){
            newn[q] = *n;
        }
        //fprintf(stdout,"calcnewN = %d\n",newn[q]);
        positioncoeff[q] = malloc((newn[q] + 1)*sizeof(double*));
        radcoeff[q] = malloc((newn[q] + 1)*sizeof(double*));
        for(int i = 0; i < (newn[q] + 1); i++){
            positioncoeff[q][i] = malloc(*dim*sizeof(double));
            radcoeff[q][i] = malloc((1)*sizeof(double));
        }

        getSplineCoeff(&positioncoeff[q],&radcoeff[q],&node0[q],&node1[q],dim,&newn[q],&findpoints[q],&comboindex[q],tolerance);
        //after each calculation we look at gradients to see if we need to update with a new one    
        //firstly add in determination if can even have a gradient :)
        //then add gradient calculations :)
        //if(hasGradient[nodeconnections[q][0]] == 0){
        //    //unclassified for this node, thus we will itterate and define with a gradient
        //    int countcc = 0;
        //    for(int i = 0; i < combocount; i++){
        //        if(nodeconnections[i][0] == nodeconnections[q][0]){
        //            //id 1 matches current look so add
        //            countcc++;
        //        }
        //        else if(nodeconnections[i][1] == nodeconnections[q][0]){
        //            //id 2 matches current look so add
        //            countcc++;
        //        }
        //    }
        //    if(countcc == 2){
        //        hasGradient[nodeconnections[q][0]] = 1;
        //        //we have marked it has a gradient as thus need to add in a gradient :)
        //        for(int i = 0; i < *dim + 1; i++){
        //            //note since nodeconnections[q][0]; we get from p1 to p0 
        //            if(i < *dim){
        //                markGradient[nodeconnections[q][0]][i] = positioncoeff[q][0][i] - positioncoeff[q][1][i];
        //            }
        //            else{
        //                markGradient[nodeconnections[q][0]][i] = radcoeff[q][0][i-*dim] - radcoeff[q][1][i-*dim];
        //            }
        //        }
        //    }
        //    else{
        //        hasGradient[nodeconnections[q][0]] = -1;
        //    }
        //}
        //if(hasGradient[nodeconnections[q][1]] == 0){
        //    //unclassified for this node, thus we will itterate and define with a gradient
        //    int countcc = 0;
        //    for(int i = 0; i < combocount; i++){
        //        if(nodeconnections[i][0] == nodeconnections[q][1]){
        //            //id 1 matches current look so add
        //            countcc++;
        //        }
        //        else if(nodeconnections[i][1] == nodeconnections[q][1]){
        //            //id 2 matches current look so add
        //            countcc++;
        //        }
        //    }
        //    if(countcc == 2){
        //        hasGradient[nodeconnections[q][1]] = 1;
        //        //we have marked it has a gradient as thus need to add in a gradient :)
        //        for(int i = 0; i < *dim + 1; i++){
        //            //note since nodeconnections[q][1]; we get from pn-1 to pn 
        //            if(i < *dim){
        //                markGradient[nodeconnections[q][1]][i] = positioncoeff[q][newn[q]][i] - positioncoeff[q][newn[q] - 1][i];
        //            }
        //            else{
        //                markGradient[nodeconnections[q][1]][i] = radcoeff[q][newn[q]][i-*dim] - radcoeff[q][newn[q] - 1][i-*dim];
        //            }
        //        }
        //    }
        //    else{
        //        hasGradient[nodeconnections[q][1]] = -1;
        //    }
        //}
        //free(applyGradient);
    }
    for(int q = 0; q < combocount; q++){
        char indxname[80];
        sprintf (indxname, "dat/splineBranchDat-%5.3f-P%d.dat", t, pid());
        FILE * fpindx = fopen (indxname, "a");
        //Out puts node0(x,y) , node1(x,y) , and then if n coeff
        if(newn[q] > 0){
            fprintf(fpindx,"%d ",newn[q]);
            for(int i = 0; i < newn[q] + 1; i++){
                fprintf(fpindx,"%f %f %f ",positioncoeff[q][i][0],positioncoeff[q][i][1],radcoeff[q][i][0]);
            }
            fprintf(fpindx,"\n");
            //fprintf(fpindx,"%f %f %f %f %f %f %f %f\n",node0[q][0],node0[q][1],node0[q][2],node1[q][0],node1[q][1],node1[q][2],positioncoeff[q][0][0],positioncoeff[q][0][1]);
        }
        fflush(fpindx);
        fclose(fpindx);
    }

    //Finally free up needed variables 
    freeCombo(&findpoints,&comboindex,&nodeconnections,&nodeindex,&combocount,&nicount);
    free(visitcount);
    for(int i = 0; i < combocount; i++){
        for(int j = 0; j < newn[i] + 1; j++){
            free(positioncoeff[i][j]);
            free(radcoeff[i][j]);
        }
        free(positioncoeff[i]);
        free(radcoeff[i]);
    }
    free(positioncoeff);
    free(radcoeff);
    free(newn);
    for(int q = 0; q < combocount; q++){
        free(node0[q]);
        free(node1[q]);
    }
    free(node0);
    free(node1);
    for(int i = 0; i < ncount; i++){
        free(nodeid[i]);
    }   
    free(nodeid);
    for(int i = 0; i < stackcount * stackcount; i++){
        free(visitstack[i]);
    }
    free(visitstack);
    for(int i = 0; i < idcount; i++){
        free(ids[i]);
    }
    free(ids);
    //for(int q = 0; q  < localcount; q++){
    //    free(markGradient[q]);
    //}
    //free(markGradient);
    //free(hasGradient);
    *sD = sDmain;
}
void skeleReduce(double **skeleton,double delta,double *minblen,int *length,int *dim,int *mxpt,double t,int n){
    if(*length > 0){
        //We will have to brute force our data, however we will target area of data max & mins
        double xmax = -HUGE;
        double xmin = HUGE;
        double ymax = -HUGE;
        double ymin = HUGE;
        double zmax = -HUGE;
        double zmin = HUGE;
        for(int i = 0; i < *length; i++){
            if(skeleton[i][0] < xmin){
                xmin = skeleton[i][0];
            }
            if(skeleton[i][0] > xmax){
                xmax = skeleton[i][0];
            }
            if(skeleton[i][1] < ymin){
                ymin = skeleton[i][1];
            }
            if(skeleton[i][1] > ymax){
                ymax = skeleton[i][1];
            }
            if(*dim == 3){
                if(skeleton[i][2] < zmin){
                    zmin = skeleton[i][2];
                }
                if(skeleton[i][2] > zmax){
                    zmax = skeleton[i][2];
                }
            }
        }
        xmax = xmax + delta * 1.5;
        xmin = xmin - delta * 1.5;
        ymax = ymax + delta * 1.5;
        ymin = ymin - delta * 1.5;
        if(*dim == 3){
            zmax = zmax + delta * 1.5;
            zmin = zmin - delta * 1.5;
        }
        double tolerance = 1e-5;
        //Next we create our Structures
        
        struct skeleDensity *sD;
        createSD(&sD,skeleton,length,dim,delta,xmax,xmin,ymax,ymin,zmax,zmin);
        makeNodePoint(&sD,dim);
        makeSpline(&sD,dim,&tolerance,length,t,&n);
        ///////////////////////////////////////////////////////////////////////Clearance
        //adapt_wavelet({hsd,hnpt,hlpt,hrpt},(double[]) {tolerance,tolerance,tolerance,tolerance}, sd.level,sd.level);
        //adapt_wavelet2({hsd,hnpt,hlpt,hrpt},(double[]) {tolerance,tolerance,tolerance,tolerance},calclevel,calclevel);
        //sd.sd = hsd; 
        //sd.npt = hnpt;
        //sd.lpt = hlpt;
        //sd.rpt = hrpt;
        //int tcount = 0; 
        //foreach_level_or_leaf(calclevel){
        //    sd.sd[] = 0.0;
        //    sd.npt[] = 0;
        //    if(x + Delta / 2 > sd.xmin && x - Delta / 2 < sd.xmax){
        //        if(y + Delta / 2 > sd.ymin && y - Delta / 2 < sd.ymax){
        //            //Here we are inside our relative bounds, now we will assign/count points 
        //            bool tcounted = false;
        //            for(int i = 0; i < *length; i++){
        //                if(skeleton[i][0] > x - Delta / 2 && skeleton[i][0] < x + Delta / 2){
        //                    if(skeleton[i][1] > y - Delta / 2 && skeleton[i][1] < y + Delta / 2){
        //                        if(!tcounted){
        //                            tcount++;
        //                            tcounted = true;
        //                        }
        //                        sd.lpt.x[0,0,sd.npt[]] = skeleton[i][0];
        //                        sd.lpt.y[0,0,sd.npt[]] = skeleton[i][1];
        //                        sd.rpt[0,0,sd.npt[]] = skeleton[i][2];
        //                        sd.npt[] = sd.npt[] + 1;
        //                    }
        //                }
        //            }
        //            //Next we calculate density for each inside
        //            double area = Delta * Delta;
        //            //fprintf(stdout,"area=%f\n",area);
        //            sd.sd[] = sd.npt[] / area;
        //        }
        //    }
        //}
        ////next we will find local maximums for density  and create root points where we think they should be
        //double **nodePoints = (double**)malloc(tcount * sizeof(double*));
        //fprintf(stdout,"\ndim=%d\n\n",*dim);
        //for(int tempi = 0; tempi < tcount; tempi++){
        //    nodePoints[tempi] = (double*)malloc((*dim + 1) * sizeof(double));
        //}
        //int npcount = 0;
        //foreach_level_or_leaf(calclevel){
        //    if (sd.npt[] > 0){
        //        //here we know we have points
        //        bool allowthrough = true;
        //        int nearcount = 0;
        //        bool skip = false;
        //        if(tcount > 9){
        //            for(int q = -1; q <= 1;q++){    
        //                for(int p = -1; p <= 1;p++){
        //                    if(sd.sd[] < sd.sd[q,p]){
        //                        allowthrough = false;
        //                    }
        //                    if(sd.npt[q,p] > 0 && q != 0 && p != 0){
        //                        nearcount++;
        //                    }
        //                }
        //            }
        //        }
        //        else{
        //            skip = true;
        //        }
        //        if((allowthrough && nearcount != 0)|| nearcount == 1 || skip){
        //            //We dont have enough to have a true local max, so we can just make a node point everywhere
        //            double ax = 0.;
        //            double ay = 0.;
        //            double ar = 0.;
        //            for(int i = 0; i < sd.npt[]; i++){
        //                ax += sd.lpt.x[0,0,i]; 
        //                ay += sd.lpt.y[0,0,i]; 
        //                ar += sd.rpt[0,0,i];
        //            }
        //            ax = ax / sd.npt[];
        //            ay = ay / sd.npt[];
        //            ar = ar / sd.npt[];
        //            nodePoints[npcount][0] = ax;
        //            nodePoints[npcount][1] = ay;
        //            nodePoints[npcount][2] = ar;
        //            npcount++;
        //        }
        //    }
        //}
        
        //Output variables
        //char npname[80];
        //sprintf (npname, "nodePoint-%5.3f.dat", t);
        //FILE * fpnp = fopen (npname, "w");
        //for(int i = 0; i < npcount; i++){
        //    fprintf(fpnp,"%f %f %f\n",nodePoints[i][0],nodePoints[i][1],nodePoints[i][2]);
        //}
        //fflush(fpnp);
        //fclose(fpnp);

        char boxname[80];
        sprintf (boxname, "dat/boxDat-%5.3f.dat", t);
        FILE * fpbox = fopen (boxname, "w");
        for(int i = 0; i < *sD->row;i++){
            for(int j = 0; j < *sD->col;j++){
                for(int k = 0; k < *sD->dep;k++){
                    //outputs --point, and ++point
                    struct skeleBounds sb = (sD->sB)[i][j][k];
                    fprintf(fpbox,"%f %f %f %f %f %f %f %f %f %d %d %d\n",*sb.x,*sb.y,*sb.x + *(sD->dx),*sb.y + *(sD->dy),sb.roi[0][0],sb.roi[0][1],sb.roi[1][0],sb.roi[1][1],*sb.density,*sD->row,*sD->col,*sb.leng > 0);
                }
            }
        }
        fflush(fpbox);
        fclose(fpbox);
        

        char nodename[80];
        sprintf (nodename, "dat/nodeDat-%5.3f.dat", t);
        FILE * fpnode = fopen (nodename, "w");
        for(int i = 0; i < *sD->row;i++){
            for(int j = 0; j < *sD->col;j++){
                for(int k = 0; k < *sD->dep;k++){
                    //outputs --point, and ++point
                    struct skeleBounds sb = (sD->sB)[i][j][k];
                    if(*sb.hasNode){
                        fprintf(fpnode,"%f %f %f %d\n",sb.nodepoint[0],sb.nodepoint[1],sb.nodepoint[2],*sb.smode);
                    }
                }
            }
        }
        fflush(fpnode);
        fclose(fpnode);
        
        char scname[80];
        sprintf (scname, "dat/splinecalcDat-%5.3f.dat", t);
        FILE * fpsc = fopen (scname, "w");
        for(int i = 0; i < *sD->row;i++){
            for(int j = 0; j < *sD->col;j++){
                for(int k = 0; k < *sD->dep;k++){
                    //outputs --point, and ++point
                    struct skeleBounds sb = (sD->sB)[i][j][k];
                    if(sb.closedis1 != NULL && *sb.closedis1 != -1){
                        fprintf(fpsc,"%f %f %f %f %d %d %d\n",*sb.x,*sb.y,*(sD->dx),*(sD->dy),*sb.closedis1,sb.closeid1[0],sb.closeid1[1]);
                    }
                }
            }
        }
        fflush(fpsc);
        fclose(fpsc);
        destroySD(&sD,dim);
        
        //cleanup var
        //for(int i = 0; i < tcount; i++){
        //    free(nodePoints[i]);
        //}
        //free(nodePoints);
        //nodePoints = NULL;

    }
}
