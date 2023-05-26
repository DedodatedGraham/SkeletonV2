#ifndef _skeleB_
#include "skeletize.h"
#include "../basiliskfunctions/adapt2.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
double** thinSkeleton(double **skeleton,int *length,double *alpha){
    fprintf(stdout,"oldL=%d\n",*length);
    //fprintf(stdout,"alpha=%f\n",*alpha);
    //fprintf(stdout,"b4:");
    //for(int i = 0; i < *length;i++){
    //    fprintf(stdout,"[%f,%f][%f,%f] , ",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3]);
    //}
    //fprintf(stdout,"\n");
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
                //skeleton[j-1] = skeleton[j];
            }
            int L = *length - 1;
            *length = L;
        }
    }
    if(addq){
        *length = *length + 1;
    }
    double **newskeleton = realloc(skeleton,(*length) * sizeof(double*));
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
    double *x;//By convention x,y is the Top left corner of the bounds;
    double *y;
    double *z;
    double **points;
    int *leng;
    double *density;
    bool *hasNode;
    double *nodepoint;
};
struct skeleDensity{
    //overall struct, this holds relevant values needed to be calculated &
    //our grid
    struct skeleBounds ***sB;//we have a 3D matrix of skele bounds refrenced by [row][col][dep]
    //NOTE: dep is only deeper than 1 if dim = 3
    int *row;
    int *col;
    int *dep;
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
                for(int q = 0; q < *tlen; q++){
                    (allocsb[i][j][k].points)[q] = malloc(*dim * sizeof(double*));
                    for(int p = 0; p < *dim; p++){
                        (allocsb[i][j][k].points)[q][p] = holdpoint[q][p];
                    }
                }
                free(holdpoint);//frees the extra memeory we dont actually want to keep from our sort
                //Finally we can calculate our density
                double *calcdensity = malloc(sizeof(double));
                *calcdensity = *tlen / (10000 * ((*dx) * (*dy) * (*dz)));//points are our 'mass' 
                allocsb[i][j][k].density = calcdensity;
                //fprintf(stdout,"grid [%f,%f,%f] / %d\n",*dx,*dy,*dz,*tlen);
                fprintf(stdout,"grid [%d,%d,%d] = %f\n\n",i,j,k,*calcdensity);
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
                    free(bounds->leng);
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
        free(tsd);
        *sd = NULL;
    }
    else{
        fprintf(stdout,"error NULL\n");
    }
}
void makeNodePoint(struct skeleDensity **sd){
    
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
    double tolerance = 1e-5;
    //Next we create our Structures
    
    struct skeleDensity *sD;
    createSD(&sD,skeleton,length,dim,delta,xmax,xmin,ymax,ymin,zmax,zmin);
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
                fprintf(fpbox,"%f %f %f %f\n",*sb.x,*sb.y,*sb.x + *(sD->dx),*sb.y + *(sD->dy));
        }
    }
    }
    fflush(fpbox);
    fclose(fpbox);
    destroySD(&sD,dim);
    
    //cleanup var
    //for(int i = 0; i < tcount; i++){
    //    free(nodePoints[i]);
    //}
    //free(nodePoints);
    //nodePoints = NULL;
}
#endif
