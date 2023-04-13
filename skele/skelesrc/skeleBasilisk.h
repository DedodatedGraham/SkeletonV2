#ifndef _skeleB_
#include "skeletize.h"
#include "../basiliskfunctions/adapt2.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
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
void getCoeffGE(int row, int col, double **a, double *x){
    int i,j,k;
    for(i = 0; i < row - 1; i++){
        //swap
        for(k = i + 1; k < row; k++){
            if(fabs(a[i][i])<fabs(a[k][i])){
                for(j = 0; j < col; j++){
                    double temp;
                    temp = a[i][j];
                    a[i][j] = a[k][j];
                    a[k][j] = temp;
                }
            }
        }
        //Gauss
        for(k = i + 1; k <row;k++){
            double term = a[k][i] / a[i][i];
            for(j = 0; j < col; j++){
                a[k][j] = a[k][j] - term*a[i][j];
            }
        }
    }
    //sub
    for(i = row - 1; i >= 0; i--){
        x[i]=a[i][col-1];
        for(j  = i + 1; j < col - 1; j++){
            x[i] = x[i] - a[i][j]*x[j];
        }
        x[i] = x[i]/a[i][i];
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
    double **arr = (double**)malloc((nrp+1)*sizeof(double*));
    for(int k = 0; k < nrp+1; k++){
        arr[k] = (double*)malloc(nc*sizeof(double));
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
    double **arr = (double**)malloc((nrp+1)*sizeof(double*));
    for(int k = 0; k < nrp+1; k++){
        arr[k] = (double*)malloc(nc*sizeof(double));
    }
    int grabarea = 3;//Odd number, the amount of area we want to scan around a point. Eg. if 5 we get -2,+2 from cell in each dim
    //Calculate the interface data
    int arrindx = 0;

    foreach(serial){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = (double**)malloc((grabarea * grabarea)*sizeof(double*));//allocated max amount of points
            for(int i = 0; i < grabarea*grabarea; i++){
                localSpline[i] = (double*)malloc((nc/2)*sizeof(double));
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
            double *X = (double*)malloc((2*n+1) * sizeof(double));
            double *Y = (double*)malloc((n + 1) * sizeof(double));
            double **B = (double**)malloc((n+1) * sizeof(double*));
            double *A = (double*)malloc((n + 1) * sizeof(double));
            for(int i = 0; i < n+1; i++){
                B[i] = (double*)malloc((n+2) * sizeof(double));
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
                getCoeffGE(n+1,n+2,B,A);
                //Finally we will Get our current point
                arr[arrindx][0] = x;
                arr[arrindx][1] = 0;
                double *AP = (double*)malloc(n*sizeof(double));
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
                getCoeffGE(n+1,n+2,B,A);
                //Finally we will Get our current point
                arr[arrindx][1] = y;
                arr[arrindx][0] = 0;
                double *AP = (double*)malloc(n*sizeof(double));
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

struct skeleDensity{
    scalar sd;
    vector lpt;
    scalar npt;
    scalar rpt;
    double xmax;
    double xmin;
    double ymax;
    double ymin;
    double zmax;
    double zmin;
    int level;
};

void makeNodePoint(struct skeleDensity sd,double **pointlist,int *position, int *dim){
    //NOTE already in a foreach_leaf_or_level
}
void skeleReduce(double **skeleton,double *minblen,int *length,int *dim,struct skeleDensity sd,int *mxpt,double t){
    //first assign all points of skeleton
    //We will have to brute force our data, however we will target area of data max & mins
    sd.xmax = -HUGE;
    sd.xmin = HUGE;
    sd.ymax = -HUGE;
    sd.ymin = HUGE;
    for(int i = 0; i < *length; i++){
        if(skeleton[i][0] < sd.xmin){
            sd.xmin = skeleton[i][0];
        }
        if(skeleton[i][0] > sd.xmax){
            sd.xmax = skeleton[i][0];
        }
        if(skeleton[i][1] < sd.ymin){
            sd.ymin = skeleton[i][1];
        }
        if(skeleton[i][1] > sd.ymax){
            sd.ymax = skeleton[i][1];
        }
    }
    
    scalar hsd[] = sd.sd;
    scalar hnpt[] = sd.npt;
    scalar hrpt[] = sd.rpt;
    vector hlpt[] = sd.lpt;
    double tolerance = 1e-3;
    //adapt_wavelet({hsd,hnpt,hlpt,hrpt},(double[]) {tolerance,tolerance,tolerance,tolerance}, sd.level,sd.level);
    adapt_wavelet2({hsd,hnpt,hlpt,hrpt},(double[]) {tolerance,tolerance,tolerance,tolerance}, sd.level+1,sd.level+1);
    sd.sd = hsd; 
    sd.npt = hnpt;
    sd.lpt = hlpt;
    sd.rpt = hrpt;
    int tcount = 0; 
    foreach_level_or_leaf(sd.level){
        sd.sd[] = 0.0;
        sd.npt[] = 0;
        if(x + Delta / 2 > sd.xmin && x - Delta / 2 < sd.xmax){
            if(y + Delta / 2 > sd.ymin && y - Delta / 2 < sd.ymax){
                //Here we are inside our relative bounds, now we will assign/count points 
                bool tcounted = false;
                for(int i = 0; i < *length; i++){
                    if(skeleton[i][0] > x - Delta / 2 && skeleton[i][0] < x + Delta / 2){
                        if(skeleton[i][1] > y - Delta / 2 && skeleton[i][1] < y + Delta / 2){
                            if(!tcounted){
                                tcount++;
                                tcounted = true;
                            }
                            sd.lpt.x[0,0,sd.npt[]] = skeleton[i][0];
                            sd.lpt.y[0,0,sd.npt[]] = skeleton[i][1];
                            sd.rpt[0,0,sd.npt[]] = skeleton[i][2];
                            sd.npt[] = sd.npt[] + 1;
                        }
                    }
                }
                //Next we calculate density for each inside
                double area = Delta * Delta;
                //fprintf(stdout,"area=%f\n",area);
                sd.sd[] = sd.npt[] / area;
            }
        }
    }
    //next we will find local maximums for density  and create root points where we think they should be
    double **nodePoints = (double**)malloc(tcount * sizeof(double*));
    for(int tempi = 0; tempi < tcount; tempi++){
        nodePoints[tempi] = (double*)malloc((*dim + 1) * sizeof(double));
    }
    int npcount = 0;
    foreach_level_or_leaf(sd.level){
        if (sd.npt[] > 0){
            //here we know we have points
            bool allowthrough = true;
            int nearcount = 0;
            if(tcount >=9){
                for(int q = -1; q <= 1;q++){    
                    for(int p = -1; p <= 1;p++){
                        if(sd.sd[] < sd.sd[q,p]){
                            allowthrough = false;
                        }
                        if(sd.npt[q,p] > 0 && q != 0 && p != 0){
                            nearcount++;
                        }
                    }
                }
            }
            if((allowthrough && nearcount != 0)|| nearcount == 1){
                //We dont have enough to have a true local max, so we can just make a node point everywhere
                double ax = 0.;
                double ay = 0.;
                double ar = 0.;
                for(int i = 0; i < sd.npt[]; i++){
                    ax += sd.lpt.x[0,0,i]; 
                    ay += sd.lpt.y[0,0,i]; 
                    ar += sd.rpt[0,0,i];
                }
                ax = ax / sd.npt[];
                ay = ay / sd.npt[];
                ar = ar / sd.npt[];
                nodePoints[npcount][0] = ax;
                nodePoints[npcount][1] = ay;
                nodePoints[npcount][2] = ar;
                npcount++;

            
            }
        }
    }
    char npname[80];
    sprintf (npname, "nodePoint-%5.3f.dat", t);
    FILE * fpnp = fopen (npname, "w");
    for(int i = 0; i < npcount; i++){
        fprintf(fpnp,"%f %f %f\n",nodePoints[i][0],nodePoints[i][1],nodePoints[i][2]);
    }
    fflush(fpnp);
    fclose(fpnp);

    char boxname[80];
    sprintf (boxname, "boxDat-%5.3f.dat", t);
    FILE * fpbox = fopen (boxname, "w");
    foreach_level_or_leaf(sd.level){
        if(x + Delta / 2 > sd.xmin && x - Delta / 2 < sd.xmax){
            if(y + Delta / 2 > sd.ymin && y - Delta / 2 < sd.ymax){
                //outputs --point, and ++point
                fprintf(fpbox,"%f %f %f %f\n",x - Delta/2, y - Delta/2, x + Delta/2, y + Delta/2);
            }
        }
    }
    fflush(fpbox);
    fclose(fpbox);

    for(int i = 0; i < tcount; i++){
        free(nodePoints[i]);
    }
    free(nodePoints);
    nodePoints = NULL;
}
#endif
