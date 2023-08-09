#ifndef _skeleB_
#include "skeletize.h"
#include "../basiliskfunctions/adapt2.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
void thinSkeleton(double ***pskeleton,int *dim,int *length,double *alpha){
    double **skeleton = *pskeleton;
    fprintf(stdout,"oldL=%d\n",*length);
    //fprintf(stdout,"alpha=%f\n",*alpha);
    //fprintf(stdout,"b4:");
    //for(int i = 0; i < *length;i++){
    //    fprintf(stdout,"[%f,%f][%f,%f] , ",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]);
    //}
    //fprintf(stdout,"\n");
    int holdl = *length;
    bool addq = false;
    for(int i = *length - 1; i >=0; i--){
        //fprintf(stdout,"alpha%d=%f\n",i,skeleton[i][3]);
        if(skeleton[i][3] < *alpha){
            //If unoptimal we will shift everything down one
            if(!addq){
                addq = true;
            }
            for(int j = i+1; j < *length; j++){
                skeleton[j-1][0] = skeleton[j][0];
                skeleton[j-1][1] = skeleton[j][1];
                skeleton[j-1][2] = skeleton[j][2];
                skeleton[j-1][3] = skeleton[j][3];
                if(*dim == 3){
                    skeleton[j-1][4] = skeleton[j][4];
                }
                //skeleton[j-1] = skeleton[j];
            }
            int L = *length - 1;
            *length = L;
        }
    }
    if(addq){
        *length = *length + 1;
    }
    double **newskeleton = malloc((*length) * sizeof(double*));
    for(int i = 0; i < *length; i++){
        newskeleton[i] = calloc(*dim + 2,sizeof(double));
        newskeleton[i][0] = skeleton[i][0];
        newskeleton[i][1] = skeleton[i][1];
        newskeleton[i][2] = skeleton[i][2];
        newskeleton[i][3] = skeleton[i][3];
        if(*dim == 3){
            newskeleton[i][4] = skeleton[i][4];
        }
    }
    for(int i = 0; i < holdl + 1; i++){
        free(skeleton[i]);
    }
    free(skeleton);
    *pskeleton = newskeleton;
    fprintf(stdout,"newL=%d\n",*length);
    //fprintf(stdout,"af:");
    //for(int i = 0; i < *length;i++){
    //    fprintf(stdout,"[%f,%f][%f,%f] , ",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]);
    //}
    //fprintf(stdout,"\n");
}


struct OutputXYNorm{
    scalar c;
    //FILE *fp;
    face vector s;
    int level;
};
//Method For ordering an interface, currently 2D only
//double** orderInterface(double **inputPoints,struct kdleaf *kdstruct,int *length,int *dim){
//    fprintf(stdout,"starting process\n");
//    //Inputs are unordered points which wont work for smoothing
//    double **newInterface = (double **)malloc(*length * sizeof(double*));//We attach new ordered points here
//    double **visitedStack = (double **)malloc(*length * sizeof(double*));//Holds the visisted points to send to kdtree
//    int j = 0;
//    while(j < *length){
//        newInterface[j] = (double*)malloc(*dim * sizeof(double));
//        j ++;
//    }
//    fprintf(stdout,"alloc\n");
//    
//    //loop var
//    //we will start at the first point and go on from there
//    int i = 0;//index for current position of point
//    double *currentpoint;
//    bool runningloop = true;
//    fprintf(stdout,"starting loop\n");
//    while(runningloop){
//        fprintf(stdout,"\non point:%d/%d \n",i,*length);
//        //running loop will complete when all points have been ordered in some way or decided to be forgotten if necessiary
//        if(i == 0){
//            //If first itteration then we will choose the first point, add it to stack and move on
//            currentpoint = inputPoints[i];
//            visitedStack[i] = inputPoints[i];
//            newInterface[i][0] = inputPoints[i][0];
//            newInterface[i][1] = inputPoints[i][1];
//            if(*dim == 3){
//                newInterface[i][2] = inputPoints[i][2];
//            }
//            i++;
//        }
//        else{
//            double lowestdistance = 0.;
//            int tempi = i;
//            double *temppoint = getNearest(currentpoint,kdstruct,length,dim,visitedStack,&tempi,&lowestdistance);
//            fprintf(stdout,"final dist %f\n",lowestdistance);
//            if(lowestdistance != 0.0){
//                //When we have conclusive results as there is a close point which has been uncounted
//                currentpoint = temppoint;//set new value
//                visitedStack[i] = temppoint;//assins old value into stack
//                newInterface[i][0] = currentpoint[0];
//                newInterface[i][1] = currentpoint[1];
//                if(*dim == 3){
//                    newInterface[i][2] = currentpoint[2];
//                }
//                i++;//increments as were moving on
//            }
//            else{
//                fprintf(stdout,"caution empty return\n");
//                fprintf(stdout,"unfound point:[%f,%f]\n",currentpoint[0],currentpoint[1]);
//                for(int j = 0; j < i; j++){
//                    fprintf(stdout,"ignoring point:[%f,%f]\n",visitedStack[j][0],visitedStack[j][1]);
//                }
//                for(int j = 0; j < i; j++){
//                    fprintf(stdout,"new:[%f,%f]\n",newInterface[j][0],newInterface[j][1]);
//                }
//                break;
//            }
//        }
//        if(i == *length){
//            //we have counted all the points
//            runningloop = false;
//        }
//    }
//    //in future full implementation since we malloc new list to seperate dependencies, 
//    //we will also free our input points as it wont get cleaned anywhere else
//    //free(inputPoints)
//    return newInterface;
//}
//
//void output_orderedlistxy(double **list,int *length,char path[80]){ 
//    FILE * fp1 = fopen (path, "w");
//    //output_facets (f,fp1);
//    for(int i = 0; i < *length; i++){
//        fprintf(fp1,"%g %g %d\n",list[i][0],list[i][1],i);
//    }
//    fflush(fp1);
//    fclose(fp1);
//}
//void smooth(double **originalInterface,int *length, int *dim,double t){
//    //For our smoothing function we are input the VOF points
//    //We will first create a robust spline of our data
//    //From the robust spline we will place the amount of original points all equally spaced along the spline
//    //Resulting in similar quality skeletons for given data input
//    //However there will be an ability to hard set the amount of points if a higher/lower resolution is desired
//    //first we will create a kdtree and order the points
//    
//    //first create tree
//    struct kdleaf *kdstruct;
//    CreateStructure(originalInterface,&kdstruct,0,dim,0,*length);//make kd-struct
//    fprintf(stdout,"\n");
//    fprintf(stdout,"made kdtree\n");
//    //order the points
//    double **newpoints;
//    newpoints = orderInterface(originalInterface,kdstruct,length,dim);
//    fprintf(stdout,"ordered\n");
//    char path[80];
//    sprintf(path,"orderedinterface-%5.3f",t);
//    output_orderedlistxy(newpoints,length,path);
//    //Apply Spline to points
//    //double **splinepoints;
//    
//    //clean up
//    kdDestroy(kdstruct);
//
//    //return final points
//    //return splinepoints;
//}

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
            double term = (*a)[k][i] / (*a)[i][i];
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
        (*x)[i] = (*x)[i] / (*a)[i][i];
    }
}

//outputs 2D x,y,norm data
double** output_points_xynorm(struct OutputXYNorm p, int *nrow,int *ndim){
    *ndim = 2;
    scalar c = p.c;
    restriction({c});
    face vector s = p.s;
    //if(!p.fp) p.fp = stdout;
    if(!s.x.i) s.x.i = -1;
    int j = 0;// number of interfacial cells 
    foreach(serial){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            j++;
	    }
    }
    int nr = j; int nc = 4;// nc is the number of column, we initialize it with 4 because we will stor x,y norm data in those columns
    int nrp = 2*nr - 1;
    *nrow = nrp;
    j = 0;
    double **arr = malloc((nrp+1)*sizeof(double*));
    for(int k = 0; k < nrp+1; k++){
        arr[k] = malloc(nc*sizeof(double));
    }
    //Calculate the interface data
    foreach(serial){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            coord n = facet_normal(point, c, s);
	        double alpha = plane_alpha(c[], n);
	        coord pc;
	        double area = plane_area_center(n, alpha, &pc);
	        if(area==0){
	            fprintf(stdout,"Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	        }
	        arr[j][0] = x+Delta*pc.x; 
	        arr[j][1] = y+Delta*pc.y;
	        double abs = sqrt(pow(n.x,2)+pow(n.y,2));
            double tx = n.x/abs;
            double ty = n.y/abs;
            arr[j][2] = tx;  
            arr[j][3] = ty;
	        arr[j+1][0] = (x+Delta*pc.x); 
	        arr[j+1][1] = -(y+Delta*pc.y);
            arr[j+1][2] = tx;  
            arr[j+1][3] = -ty;
	        j = j + 2;
	    }
    }
    return arr;
}

double** output_points_2smooth(struct OutputXYNorm p, int *nrow,int *ndim, double t){
    *ndim = 2;
    scalar c = p.c;
    restriction({c});
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    int pj = 0;// number of interfacial cells
    scalar vofx[],vofy[],normx[],normy[];
    double xmax = 0.,xmin = 0.,ymax = 0.,ymin = 0.;
    foreach(serial){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            coord n = facet_normal(point, c, s);
	        double alpha = plane_alpha(c[], n);
	        coord pc;
	        double area = plane_area_center(n, alpha, &pc);
	        if(area==0){
	            fprintf(stdout,"Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	        }
            vofx[] = (x+Delta*pc.x);
            vofy[] = (y+Delta*pc.y);
            normx[] = n.x;
            normy[] = n.y;
            //Here we capture the min and max's to prevent incorrect points from making it to the scheme when we use a grid larger than 3x3
            if(vofx[] > xmax || xmax == 0.){
                xmax = vofx[];
            }
            if(vofx[] < xmin || xmin == 0.){
                xmin = vofx[];
            }
            if(vofy[] > ymax || ymax == 0.){
                ymax = vofy[];
            }
            if(vofy[] < ymin || ymin == 0.){
                ymin = vofy[];
            }
            pj++;
	    }
    }
    int nr = pj; int nc = 4;// nc is the number of column, we initialize it with 4 because we will stor x,y norm data in those columns
    int nrp = nr - 1;
    *nrow = nrp;
    double **arr = malloc((nrp+1)*sizeof(double*));
    for(int k = 0; k < nrp+1; k++){
        arr[k] = malloc(nc*sizeof(double));
    }
    int grabarea = 3;//Odd number, the amount of area we want to scan around a point. Eg. if 5 we get -2,+2 from cell in each dim
    //Calculate the interface data
    int arrindx = 0;

    foreach(serial){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = malloc((grabarea * grabarea)*sizeof(double*));//allocated max amount of points
            for(int i = 0; i < grabarea*grabarea; i++){
                localSpline[i] = malloc((nc/2)*sizeof(double));
            }
            //First we will go though and collect all needed 
            //First in X
            int indx = 0;
            for(int i = -1 * (grabarea - 1) / 2 ;i <= (grabarea - 1) / 2;i++){
                //Then in Y
                for(int j = - 1 * (grabarea - 1) / 2;j <= (grabarea - 1) / 2;j++){
                    //we can now go through each stencil and find all interface points along the area
                    if(c[i,j] > 1e-6 && c[i,j] < 1.-1e-6){
                        if(vofx[i,j] <= xmax && vofx[i,j] >= xmin && vofy[i,j] <= ymax && vofy[i,j] >= ymin){
                            localSpline[indx][0] = vofx[i,j];
                            localSpline[indx][1] = vofy[i,j];
                            indx++;
                            
                        }
                    }
                }
            }
            //Now we have all the points we want, so we next find our approx valuesj
            int n = 2;
            //allocate
            double *X = (double*)calloc((2*n+1) , sizeof(double));
            double *Y = (double*)calloc((n + 1) , sizeof(double));
            double **B = (double**)calloc((n+1) , sizeof(double*));
            double *A = (double*)calloc((n + 1) , sizeof(double));
            for(int i = 0; i < n+1; i++){
                B[i] = malloc((n+2) * sizeof(double));
            }
            //calc arrays
            if(fabs(normy[]) > fabs(normx[])){ 
                for(int i = 0; i <= 2*n; i++){
                    X[i] = 0;
                    for(int j = 0; j < indx; j++){
                        X[i] = X[i] + pow(localSpline[j][0],i);
                    }
                }
                for(int i = 0; i <= n; i++){
                    Y[i] = 0;
                    for(int j = 0; j < indx; j ++){
                        Y[i] = Y[i] + pow(localSpline[j][0],i)*localSpline[j][1];
                    }
                }
                //make B
                for(int i = 0; i <= n; i++){
                    for(int j = 0; j <= n; j++){
                        B[i][j] = X[i+j];
                    }
                }
                for(int i = 0; i <= n; i++){
                    B[i][n+1] = Y[i];
                }
                getCoeffGE(n+1,n+2,&B,&A);
                //Finally we will Get our current point
                arr[arrindx][0] = x;
                arr[arrindx][1] = 0;
                double *AP = (double*)calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    arr[arrindx][1] = arr[arrindx][1] + pow(x,i) * A[i];
                    AP[i] = A[i] * i;
                }
                //and then calculate the norms using the prime
                //First we calculate the tangent m at our point
                double m = 0;
                for(int i = 0; i <= n; i++){
                    if(i != 0){
                        m = m + AP[i] * pow(x,i-1);
                    }
                }
                //normal m = -1/m
                m = -1 * (1/m);
                double b = (-1 * m * arr[arrindx][0]) + arr[arrindx][1];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double tx;
                if(normx[] < 0.){
                    tx = arr[arrindx][0] - 1;
                }
                else{
                    tx = arr[arrindx][0] + 1;
                }
                double ty = m * tx + b;
                //Find direction vector to make
                double tnormx = tx - arr[arrindx][0];
                double tnormy = ty - arr[arrindx][1];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                arr[arrindx][2] = tnormx / bottom;
                arr[arrindx][3] = tnormy / bottom;
                //ensure norms are correct direction before outputting
                if(arr[arrindx][2] > 0. && !(normx[] > 0.)){
                    arr[arrindx][2] = -1 * arr[arrindx][2];
                }
                else if(arr[arrindx][2] < 0. && !(normx[] < 0.)){
                    arr[arrindx][2] = -1 * arr[arrindx][2];
                }
                if(arr[arrindx][3] > 0. && !(normy[] > 0.)){
                    arr[arrindx][3] = -1 * arr[arrindx][3];
                }
                else if(arr[arrindx][3] < 0. && !(normy[] < 0.)){
                    arr[arrindx][3] = -1 * arr[arrindx][3];
                }
                free(A);
                free(AP);
            }
            else{
                for(int i = 0; i <= 2*n; i++){
                    X[i] = 0;
                    for(int j = 0; j < indx; j++){
                        X[i] = X[i] + pow(localSpline[j][1],i);
                    }
                }
                for(int i = 0; i <= n; i++){
                    Y[i] = 0;
                    for(int j = 0; j < indx; j ++){
                        Y[i] = Y[i] + pow(localSpline[j][1],i)*localSpline[j][0];
                    }
                }
                //make B
                for(int i = 0; i <= n; i++){
                    for(int j = 0; j <= n; j++){
                        B[i][j] = X[i+j];
                    }
                }
                for(int i = 0; i <= n; i++){
                    B[i][n+1] = Y[i];
                }
                getCoeffGE(n+1,n+2,&B,&A);
                //Finally we will Get our current point
                arr[arrindx][1] = y;
                arr[arrindx][0] = 0;
                double *AP = (double*)calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    arr[arrindx][0] = arr[arrindx][0] + pow(y,i) * A[i];
                    AP[i] = A[i] * i;
                }
                //and then calculate the norms using the prime
                //First we calculate the tangent m at our point
                double m = 0;
                for(int i = 0; i <= n; i++){
                    if(i != 0){
                        m = m + AP[i] * pow(y,i-1);
                    }
                }
                //normal m = -1/m
                m = -1 * (1/m);
                double b = (-1 * m * arr[arrindx][1]) + arr[arrindx][0];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double ty;
                if(normy[] < 0.){
                    ty = arr[arrindx][1] - 1;
                }
                else{
                    ty = arr[arrindx][1] + 1;
                }
                double tx = m * ty + b;
                //Find direction vector to make
                double tnormx = tx - arr[arrindx][0];
                double tnormy = ty - arr[arrindx][1];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                arr[arrindx][2] = tnormx / bottom;
                arr[arrindx][3] = tnormy / bottom;
                //ensure normals are pointing in correct direction
                if(arr[arrindx][2] > 0. && !(normx[] > 0.)){
                    arr[arrindx][2] = -1 * arr[arrindx][2];
                }
                else if(arr[arrindx][2] < 0. && !(normx[] < 0.)){
                    arr[arrindx][2] = -1 * arr[arrindx][2];
                }
                if(arr[arrindx][3] > 0. && !(normy[] > 0.)){
                    arr[arrindx][3] = -1 * arr[arrindx][3];
                }
                else if(arr[arrindx][3] < 0. && !(normy[] < 0.)){
                    arr[arrindx][3] = -1 * arr[arrindx][3];
                }
                free(A);
                free(AP);
            }
            
            //Here we save data for each individual point
            //First will be VOF
            //[x,y]  [nx,ny]
            //then smooth
            //[sx,sy] [snx,sny]
            //in [x y nx ny sx sy snx sny]
            char tempfile[80];
            sprintf (tempfile, "intdata-%5.3f.dat", t);
            FILE *fp1 = fopen(tempfile,"a");
            double tempnormabs = sqrt(pow(normx[],2)+pow(normy[],2));
            fprintf(fp1,"%f %f %f %f %f %f %f %f\n",vofx[],vofy[],normx[]/tempnormabs,normy[]/tempnormabs,arr[arrindx][0],arr[arrindx][1],arr[arrindx][2],arr[arrindx][3]);
            fflush(fp1);
            fclose(fp1);
            arrindx = arrindx + 1;
            //freeup variables
            for(int i = 0;i < grabarea * grabarea; i++){
                free(localSpline[i]);
            }
            free(localSpline);
            for(int i = 0; i < n+1; i++){
                free(B[i]);
            }
            free(B);
            free(X);
            free(Y);
        }
    }
    return arr;
}
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
    fprintf(stdout,"complete\n");
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
        fprintf(stdout,"error NULL\n");
    }
}
bool connectingROI(double **roi1,double**roi2,int *dim,double *tolerance){
    //roi is [xmin,ymin,zmin],[xmax,ymax,zmax]
    //a connection has atleast one touching dimension
    bool ret = false;
    if(*dim == 2){
        //ignore z bc will always be touching
        //First we check if one dimension lines up, this determines if they are even relative to eachother
        //then check if other dimension is between
        if((abs(roi1[0][0]-roi2[0][0])<*tolerance)||(abs(roi1[0][0]-roi2[1][0])<*tolerance)||(abs(roi1[1][0]-roi2[0][0])<*tolerance)||(abs(roi1[1][0]-roi2[1][0])<*tolerance)){
            if((abs(roi1[0][1]-roi2[0][1])<*tolerance)||(abs(roi1[0][1]-roi2[1][1])<*tolerance)||(abs(roi1[1][1]-roi2[0][1])<*tolerance)||(abs(roi1[1][1]-roi2[1][1])<*tolerance)){
                ret = true;
            }
        }
    }
    return ret;
}
//bool passcheck(int i, int j, int k, int *dim){
//    if(*dim == 2){
//        if((i == 0 && j != 0) || (i != 0 && j == 0)){
//            return true;
//        }
//    }
//    else{
//        if((i == 0 && j == 0 && k != 0) || (i == 0 && j != 0 && k == 0) || (i != 0 && j == 0 && k == 0)){
//            return true;
//        }
//    }
//    return false;
//}
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
                    int nodecount = 0;
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
                        double cz;
                        if(*dim == 3){
                            cz = (*(sDmain->zmin) + *(sDmain->zmax)) / 2;
                        }
                        //And we will select the furthest point from the 'center'
                        double furthest = 0.;
                        int holdq = 0;
                        for(int q = 0; q < *(localcell->leng); q++){
                            double dx = localcell->points[q][0] - cx;
                            double dy = localcell->points[q][1] - cy;
                            double d;
                            if(*dim == 3){
                                double dz = localcell->points[q][2] - cz;
                                d = pow(dx,2) + pow(dy,2) + pow(dz,2); 
                            }
                            else{
                                d = pow(dx,2) + pow(dy,2); 
                            }
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
                        double az;
                        if(*dim == 3){
                            az = 0.;
                        }
                        //fprintf(stdout,"leng=%d\n",*(localcell->leng));
                        for(int q = 0; q < *(localcell->leng); q++){
                            ax += localcell->points[q][0];
                            ay += localcell->points[q][1];
                            if(*dim == 3){
                                az += localcell->points[q][2];
                                ar += localcell->points[q][3];
                            }
                            else{
                                ar += localcell->points[q][2];
                            }
                            //fprintf(stdout,"summing[%f,%f,%f]\n",ax,ay,ar);
                        }
                        ax = ax / *(localcell->leng);
                        ay = ay / *(localcell->leng);
                        if(*dim == 3){
                            az = az / *(localcell->leng);
                        }
                        ar = ar / *(localcell->leng);
                        calcpoint[0] = ax;
                        calcpoint[1] = ay;
                        if(*dim == 3){
                            calcpoint[2] = az;
                            calcpoint[3] = ar;
                        }
                        else{
                            calcpoint[2] = ar;
                        }
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
        nialloc++;
    }
    if(ic == -1){
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
        fprintf(stdout,"comboadd error !\n");
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
                                        double nangle;
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
sprintf (conname, "connectionDat-%5.3f.dat", t);
    FILE * fpcon = fopen (conname, "w");
    for(int i = 0; i < *combocount;i++){
        fprintf(fpcon,"%d %d\n",nodeconnections[i][0],nodeconnections[i][1]);
    }
    fflush(fpcon);
    fclose(fpcon);
    //Connection Index:
    char indxname[80];
    sprintf (indxname, "connectionidDat-%5.3f.dat", t);
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
        double *nearpoint = getNearest(comppoints[i],&searchstruct,&lengpoints,dim,&ignorepoint,&tleng,&lowestdis);
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
                CreateStructure(ntpoints,&kdstruct_petp,0,dim,0,splinecalc,1);
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
                CreateStructure(ntpoints,&kdstruct_petn,0,dim,0,splinecalc,1);
                error_neg = calcBezierErr(kdstruct_petn,findpoints,comboindex,dim);
                kdDestroy(&kdstruct_petn);
                //free points
                for(int k = 0; k < splinediv; k++){
                    free(ntpoints[k]);
                }
                free(ntpoints);
                double der;
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
        fprintf(stdout,"itt:%d|l:%f|e:%f|el:%f|ell:%f \n",itt,lambda,error,error_last,error_lastlast);
        error_lastlast = error_last;
        error_last = error;
        itt++;
    }
    fprintf(stdout,"finished at:%d\n",itt);
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
//factorial
double getFactorial(int num){
    int ret = 1;
    for(int i = 2; i < num; i++){
        ret *= ret;
    }
    return (double)ret;
}
//Least square fitting splines
void findBestFit2(double ***ppositioncoeff,double ***pradcoeff,double **findpoints,int comboindex,int *dim, int *n,double *tolerance){
    //note can only handle cubic splines as of now
    double **positioncoeff = *ppositioncoeff;
    double **radcoeff = *pradcoeff;
    //using least square to solve for spline
    double *B = calloc(*n + 1,sizeof(double));
    //collection terms
    double **coldims = malloc((*n + 1) * sizeof(double*));
    for(int i = 0; i < *n + 1; i++){
        coldims[i] = calloc(*dim + 1,sizeof(double));
    }
    double **A = malloc(*n * sizeof(double*));
    for(int i = 0; i < *n; i++){
        A[i] = calloc(*n, sizeof(double));
    }
    //construct t for finding values
    double *t = malloc(comboindex * sizeof(double));
    for(int i = 0; i < comboindex; i++){
        //for each data point try and guess the t value of the data based on linear sampling
        //allows for least square method to match up best guess spline points to data better based on flucuations
        double tx = (positioncoeff[0][0] - findpoints[i][0]) / (positioncoeff[0][0] - positioncoeff[*n][0]);
        double ty = (positioncoeff[0][1] - findpoints[i][1]) / (positioncoeff[0][1] - positioncoeff[*n][1]);
        double tr = (radcoeff[0][0] - findpoints[i][2]) / (radcoeff[0][0] - radcoeff[*n][0]);
        if(*dim == 2){
            t[i] = (tx + ty + tr) / 3;
        }
        else{
            fprintf(stdout,"error not implemented-3D :(((( \n");
        }
        //t[i] = (double)i / (comboindex - 1);
    }
    //Next we calc using least square partials :)
    for(int i = 0; i < comboindex; i++){
        double tt = 1. - t[i];
        //Find B values at current points approx t
        for(int j = 0; j < (*n + 1); j++){
            B[j] = pow(tt,*n) * (t[i],*n) * (getFactorial(*n) / (getFactorial(j) * getFactorial(*n - j)));
        }
        for(int j = 1; j < *n; j++){
            for(int q = 0; q < *dim; q++){
                coldims[j][q] += B[j] * (findpoints[i][q] - B[0] * positioncoeff[0][q] - B[*n] * positioncoeff[*n][q]);
            }
            coldims[j][*dim] += B[j] * (findpoints[i][*dim] - B[0] * radcoeff[0][0] - B[*n] * radcoeff[*n][0]);
        }
        for(int j = 0; j < *n + 1; j++){
            for(int q = 0; q < *n + 1; q++){
                //Sum A terms for common denom
                A[j][q] += B[j] * B[q];
            }
        }
    }
    //Next we take the determinate of A
    double denom = 0.;

    //finally solve for spline points
    for(int i = 1; i < *n; i++){

    }
    //for(int i = 0; i < *dim; i++){
    //    positioncoeff[1][i] = (A2*C1[i] - A12 * C2[i]) / denom;
    //    positioncoeff[2][i] = (A1*C2[i] - A12 * C1[i]) / denom;
    //}
    //radcoeff[1][0] = (A2*C1[*dim] - A12 * C2[*dim]) / denom;
    //radcoeff[2][0] = (A1*C2[*dim] - A12 * C1[*dim]) / denom;
    //Free
    free(t);
    for(int i = 0; i < (*n + 1); i++){
        free(C[i]);
    }
    free(C);
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
            findBestFit2(&positioncoeff,&radcoeff,findpoints,comboindex,dim,n,tolerance);
        }
        else{
            fprintf(stdout,"error too little input points defaulting spline to linear/input\n");
            if(comboindex == 1){
                for(int i = 1; i < *n; i++){
                    positioncoeff[i][0] = findpoints[0][0];
                    positioncoeff[i][1] = findpoints[0][1];
                    radcoeff[i][0] = findpoints[0][2];
                }
            }
            else{
                int mid = (int)floor((*n + 1) / 2);
                for(int i = 1; i < *n; i++){
                    if(i < mid){
                        positioncoeff[i][0] = findpoints[0][0];
                        positioncoeff[i][1] = findpoints[0][1];
                        radcoeff[i][0] = findpoints[0][2];
                    }
                    else{
                        positioncoeff[i][0] = findpoints[comboindex - 1][0];
                        positioncoeff[i][1] = findpoints[comboindex - 1][1];
                        radcoeff[i][0] = findpoints[comboindex - 1][2];
                    }
                }

            }
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
    fprintf(stdout,"nodepoints b4 = %d\n",ncount);
    if(*(sDmain->row) > 7 || *(sDmain->col) > 7 || *(sDmain->dep) > 7){
        reduceLocalMax(sD,dim);
    }
    fprintf(stdout,"nodepoints af = %d\n",*sDmain->ncount);
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
    //Now weve assigned nodes & collected splines approx points
    //So we will go through and calculate our approximations
    //first we define some things we want
    //given we want order 'n' we allocate accordingly
    
    //allocation for coeffs
    double ***positioncoeff = malloc(combocount*sizeof(double**));
    double ***radcoeff = malloc(combocount*sizeof(double**));//NOTE rad coeff is set to be build like position coeff to allow approximation of other var
    double **node0 = malloc(combocount*sizeof(double*));
    double **node1 = malloc(combocount*sizeof(double*));
    for(int q = 0; q < combocount; q++){
        positioncoeff[q] = malloc((*n+1)*sizeof(double*));
        radcoeff[q] = malloc((*n+1)*sizeof(double*));
        for(int i = 0; i < (*n + 1); i++){
            positioncoeff[q][i] = malloc(*dim*sizeof(double));
            radcoeff[q][i] = malloc((1)*sizeof(double));
        }
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
        getSplineCoeff(&positioncoeff[q],&radcoeff[q],&node0[q],&node1[q],dim,n,&findpoints[q],&comboindex[q],tolerance);
    }
    for(int q = 0; q < combocount; q++){
        char indxname[80];
        sprintf (indxname, "splineBranchDat-%5.3f.dat", t);
        FILE * fpindx = fopen (indxname, "a");
        //Out puts node0(x,y) , node1(x,y) , and then if n coeff
        if(*n == 1){
            fprintf(fpindx,"%f %f %f %f %f %f %f %f\n",node0[q][0],node0[q][1],node0[q][2],node1[q][0],node1[q][1],node1[q][2],positioncoeff[q][0][0],radcoeff[q][0][0]);
        }
        else if(*n > 1){
            fprintf(fpindx,"%d ",*n);
            for(int i = 0; i < *n + 1; i++){
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
        for(int j = 0; j < *n+1; j++){
            free(positioncoeff[i][j]);
            free(radcoeff[i][j]);
        }
        free(positioncoeff[i]);
        free(radcoeff[i]);
    }
    free(positioncoeff);
    free(radcoeff);
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
    *sD = sDmain;
}
void skeleReduce(double **skeleton,double delta,double *minblen,int *length,int *dim,int *mxpt,double t,int n){
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
    sprintf (boxname, "boxDat-%5.3f.dat", t);
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
    sprintf (nodename, "nodeDat-%5.3f.dat", t);
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
    sprintf (scname, "splinecalcDat-%5.3f.dat", t);
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
#endif
