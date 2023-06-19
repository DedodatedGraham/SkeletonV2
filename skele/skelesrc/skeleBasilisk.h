#ifndef _skeleB_
#include "skeletize.h"
#include "../basiliskfunctions/adapt2.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
double** thinSkeleton(double **skeleton,int *dim,int *length,double *alpha){
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
    skeleton = newskeleton;
    fprintf(stdout,"newL=%d\n",*length);
    //fprintf(stdout,"af:");
    //for(int i = 0; i < *length;i++){
    //    fprintf(stdout,"[%f,%f][%f,%f] , ",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]);
    //}
    //fprintf(stdout,"\n");
    return skeleton;
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
    int *closedis;//closest distance to node point
    int *closeid;//closest node's id
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
                int *tcd = malloc(sizeof(int));
                *tcd = -1;
                allocsb[i][j][k].closedis = tcd;
                int *tci = calloc(3,sizeof(int));
                allocsb[i][j][k].closeid = tci;
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
                double extra = 1;
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
                    free(bounds->closedis);
                    free(bounds->closeid);
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
bool connectingROI(double **roi1,double**roi2,int *dim){
    //roi is [xmin,ymin,zmin],[xmax,ymax,zmax]
    //a connection has atleast one touching dimension
    bool ret = false;
    for(int i = 0; i < *dim; i++){
        if(roi1[0][i] == roi2[0][i] || roi1[0][i] == roi2[1][i] || roi1[1][i] == roi2[0][i] || roi1[1][i] == roi2[1][i]){
            ret = true;
            break;
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
int floodDensity(struct skeleDensity **sD,int *dim,int *nodeLocation,int distance,int **visitStack,int *visitcount){
    //Firstly we set our curent cell's distance & assign it a new 
    struct skeleDensity *sd = *sD;
    if(distance == 0){
        visitStack[0][0] = nodeLocation[0];
        visitStack[0][1] = nodeLocation[1];
        visitStack[0][2] = nodeLocation[2];
    }
    int i = visitStack[*visitcount][0];
    int j = visitStack[*visitcount][1];
    int k = visitStack[*visitcount][2];
    //fprintf(stdout,"@[%d,%d,%d];distance=%d,vc=%d\n",i,j,k,distance,*visitcount);
    *visitcount = *visitcount + 1;
    if(*(sd->sB[i][j][k]).closedis == -1 || *(sd->sB[i][j][k]).closedis > distance){
        *(sd->sB[i][j][k]).closedis = distance;
        (sd->sB[i][j][k]).closeid[0] = nodeLocation[0];
        (sd->sB[i][j][k]).closeid[1] = nodeLocation[1];
        (sd->sB[i][j][k]).closeid[2] = nodeLocation[2];
        //fprintf(stdout,"seen node[%d,%d,%d]\n",(sd->sB[i][j][k]).closeid[0],(sd->sB[i][j][k]).closeid[1],(sd->sB[i][j][k]).closeid[2]);
    }
    else{
        return 0;
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
                                    //if(connectingROI((sd->sB[i][j][k]).roi,(sd->sB[ti][tj][tk]).roi,dim)){
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
    //this will let us know if we need to reach to corners or not
    //fprintf(stdout,"\n");
    for(int q = 0; q < len; q++){
        //fprintf(stdout,"trying for [%d,%d,%d] (%d/%d)\n",i+searchid[q][0],j+searchid[q][1],k+searchid[q][2],q+1,len);
        if((searchid[q][0] != 0 && searchid[q][1] == 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] != 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] == 0 && searchid[q][2] != 0)){
            bool allowpass = true;
            for(int p = 0; p < *visitcount; p++){
                //check if node is already visited
                //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
                if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
                    allowpass = false;
                    break;
                }
            }
            if(allowpass){
                //Here we have determined that we are next to the cell & it hasnt been visited yet, so we will send a flood
                //fprintf(stdout,"going to%d++\n",*visitcount);
                visitStack[*visitcount][0] = i + searchid[q][0];
                visitStack[*visitcount][1] = j + searchid[q][1];
                visitStack[*visitcount][2] = k + searchid[q][2];
                //*visitcount = *visitcount + 1;
                edgesvisited = edgesvisited + floodDensity(&sd,dim,nodeLocation,distance + 1,visitStack,visitcount);
            }
        }
        else{
            //fprintf(stdout,"wascorner[%d,%d,%d]\n",searchid[q][0],searchid[q][1],searchid[q][2]);
            bool allowpass = true;
            for(int p = 0; p < *visitcount; p++){
                //check if node is already visited
                //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
                if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
                    allowpass = false;
                    break;
                }
            }
            if(allowpass){
                corners++;
            }
        }
    }
    //now we correct to corners if no other options
    if(corners > 0){
        for(int q = 0; q < len; q++){
            //fprintf(stdout,"trying for [%d,%d,%d] (%d/%d)\n",i+searchid[q][0],j+searchid[q][1],k+searchid[q][2],q+1,len);
            if((searchid[q][0] != 0 && searchid[q][1] != 0 && searchid[q][2] == 0) || (searchid[q][0] == 0 && searchid[q][1] != 0 && searchid[q][2] != 0) || (searchid[q][0] != 0 && searchid[q][1] == 0 && searchid[q][2] != 0)){
                bool allowpass = true;
                for(int p = 0; p < *visitcount; p++){
                    //check if node is already visited
                    //fprintf(stdout,"visted %d-[%d,%d,%d]\n",p,visitStack[p][0],visitStack[p][1],visitStack[p][2]);
                    if((i+searchid[q][0] == visitStack[p][0]) && (j+searchid[q][1] == visitStack[p][1]) && (k+searchid[q][2] == visitStack[p][2])){
                        allowpass = false;
                        break;
                    }
                }
                if(allowpass){
                    //Here we have determined that we are next to the cell & it hasnt been visited yet, so we will send a flood
                    //fprintf(stdout,"going to%d++\n",*visitcount);
                    visitStack[*visitcount][0] = i + searchid[q][0];
                    visitStack[*visitcount][1] = j + searchid[q][1];
                    visitStack[*visitcount][2] = k + searchid[q][2];
                    //*visitcount = *visitcount + 1;
                    floodDensity(&sd,dim,nodeLocation,distance + 1,visitStack,visitcount);
                }
            }
        }
    }
    //fprintf(stdout,"\n");
    for(int q = 0; q < 26; q++){
        free(searchid[q]);
    }
    free(searchid);
    *sD = sd;
    return 1;
}
//resets all visited stacks to 0
void resetFloodStack(int **vs, int leng){
    for(int i = 0; i < leng; i++){
        vs[i][0] = 0;
        vs[i][1] = 0;
        vs[i][2] = 0;
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
                                                if(!(ic == 0 && jc == 0 && kc == 0)){
                                                    //we now have ti,tj,and tk which fall into the allowable regions(ie not off grid)
                                                    //Next we shall count the amount of nearby nodes made in smode = 1
                                                    if(*(sd->sB[ti][tj][tk]).smode == 1){
                                                        //because too close to endpoint, we disable the nodepoint
                                                        *(sd->sB[ti][tj][tk]).hasNode = false;
                                                        free((sd->sB[ti][tj][tk]).nodepoint);
                                                        *(sd->sB[ti][tj][tk]).smode = 0;
                                                    }
                                                }
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
//Splining method, for merging unessicary nodes, and creating a spline
void makeSpline(struct skeleDensity **sD,int *dim){
    //Here the only node points are Skeleton points ends, we move from these
    //along the region of interest 
    struct skeleDensity *sDmain = *sD;
    //for holding our locations & amounts of endpoints
    int ncount = *sDmain->ncount;
    int pcount = *sDmain->pcount; 
    int **nodeid = malloc(ncount * sizeof(int*));
    for(int i = 0; i < ncount; i++){
        nodeid[i] = calloc(3,sizeof(int));
    }
    //first we collect the locations of each node
    int localcount = 0;
    int stackcount = 0;//allocates max needed space in our stack
    //make adjustments if needed
    if(*(sDmain->row) > 7 || *(sDmain->col) > 7 || *(sDmain->dep) > 7){
        reduceLocalMax(sD,dim);
    }
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
    int **visitstack = malloc(stackcount * sizeof(int*));
    for(int i = 0; i < stackcount; i++){
        visitstack[i] = calloc(3,sizeof(int));
    }
    int *visitcount = calloc(1,sizeof(int));
    for(int i = 0; i < localcount; i++){
        *visitcount = 0;
        floodDensity(sD,dim,nodeid[i],0,visitstack,visitcount);
        if(i != localcount-1){
            resetFloodStack(visitstack,stackcount);
        }
    }
    //We have flooded our tree with endpoint calculations,
    //Next we want to merge all of our localpoints which have been flagged, and turn the close ones into one node point
    //now we free variables 
    free(visitcount);
    for(int i = 0; i < ncount; i++){
        free(nodeid[i]);
    }   
    free(nodeid);
    for(int i = 0; i < stackcount; i++){
        free(visitstack[i]);
    }
    free(visitstack);
}
void skeleReduce(double **skeleton,double delta,double *minblen,int *length,int *dim,int *mxpt,double t){
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
    makeSpline(&sD,dim);
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
                if(sb.closedis != NULL && *sb.closedis != -1){
                    fprintf(fpsc,"%f %f %f %f %d %d %d\n",*sb.x,*sb.y,*(sD->dx),*(sD->dy),*sb.closedis,sb.closeid[0],sb.closeid[1]);
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
