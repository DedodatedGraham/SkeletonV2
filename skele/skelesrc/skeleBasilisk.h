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
	            printf("Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
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
//to keep x,y,nx,and ny data together we create 
struct smooths{
    double **points;//will be [[x0,y0,nx0,ny0],...,[xl-1,yl-1,nxl-1,nyl-1]]
    double **smoothpoints;//will be [[x0,y0,nx0,ny0],...,[xl-1,yl-1,nxl-1,nyl-1]]
#if _MPI
    int *mpicomputed;
#endif
    //ensures easy looping through all points
    int *length;
};
attribute{
    struct smooths *smooth;
}

//functions for verifying information/sorting
int ensureDistance(double *p1,double *p2, double distance){
    double neededdis = distance * 3.;//3 * delta should be max distance of a stencils reach
#if dimension == 2
    int dim = 2;
#else
    int dim = 3;
#endif
    if(getDistance(p1,p2,&dim) > neededdis){
        return 0;
    } 
    return 1;
}
#if dimension == 2
int insidebounds(double *ppoint,double tx, double ty, double tdelta){
#else
int insidebounds(double *ppoint,double tx, double ty, double tz, double tdelta){
#endif
    double half = tdelta/2;
    double maxx = tx + half;
    double minx = tx - half;
    double maxy = ty + half;
    double miny = ty - half;
#if dimension == 3
    double maxz = tz + half;
    double minz = tz - half;
#endif
    if(ppoint[0] < maxx && ppoint[0] > minx){
        if(ppoint[1] < maxy && ppoint[1] > miny){
#if dimension == 3
            if(ppoint[2] < maxz && ppoint[2] > minz)
#endif
                return 1;
        }
    }
    return 0;
}

#if TREE
//when prolongating we use our structure to find the location of where our new id should go
static void smoothvof_prolongation(Point point, scalar s){
    double val = s[];
    if(s[] != nodata){
        double *ppoint = s.smooth->points[(int)s[]];
        foreach_child(){
            //We need to refrence the structure to determine which cell will have our refrence, otherwise we set value to -1
#if dimension == 2
            if(insidebounds(ppoint,x,y,Delta)){
#else
            if(insidebounds(ppoint,x,y,z,Delta)){
#endif
                s[] = val;
            }
            else{
                s[] = nodata;
            }
        }
    }
    else{
        foreach_child(){
            s[] = nodata;
        }
    }
}

static inline void smoothvof_restriction(Point point, scalar s){
    int pass = 0;
    double val = nodata;
    foreach_child(){//we pick the first child to hold a point when we upscale, allows for even distribution of points/smoothening
        if(!pass && s[] != nodata){
            val = s[];
            pass++;
        }
    }
    s[] = val;
}
#endif

//Basilisk -> Skeleton -> Basilisk MPI comunication 
#if _MPI
void smooth_interface_MPI(struct OutputXYNorm p,scalar vofref,double t,int max_level){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#if dimension == 2
    int nc = 4;
#else
    int nc = 6;
#endif
    scalar c = p.c;
    struct smooths *smooth = malloc(sizeof(struct smooths));
    //next we calc our total id across all MPI, to ensure we build our ID's correctly
    int *localcalc = calloc(comm_size,sizeof(int));//counts our local values and communicates them
    //this will make life vastly easier when smoothening across boundaries
    double mdel = HUGE;
    foreach(reduction(+:localcalc[:comm_size]) reduction(min:mdel)){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            localcalc[pid()]++;
            if(Delta < mdel){
                mdel = Delta;
            }
        }
    }
    int *calcl = calloc(1,sizeof(int));
    int lstart = 0;
    for(int i = 0; i < comm_size; i++){
        //printf("size %d => %d\n",i,localcalc[i]);
        if(i <= pid() - 1){
            lstart += localcalc[i];
        }
        *calcl += localcalc[i];
    }
    //printf("totalsize = %d, local size stat = %d,local size end = %d\n",*calcl,lstart,lstart + localcalc[pid()]);
    smooth->length = calcl;
    double **vofpoints = malloc(*calcl * sizeof(double*));
    double **smoothpoints = malloc(*calcl * sizeof(double*));
    int *beencalc = calloc(*calcl , sizeof(int));
    for(int i = 0; i < *calcl; i++){
        vofpoints[i] = calloc(nc,sizeof(double));
        smoothpoints[i] = calloc(nc,sizeof(double));
        beencalc[i] = -1; 
    }
    smooth->points = vofpoints;
    smooth->smoothpoints = smoothpoints;
    smooth->mpicomputed = beencalc;
    vofref.smooth = smooth;//sets the field to point at our structure
    //so now we have allocated a structure belonging to our ref which has enough size to hold all of the points, regardless of if they cross a boundary
    //however we will only assign in points which fall across our foreach_neighbor; first we allocate our ID's
    foreach(noauto){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            vofref[] = lstart++; 
        }
        else{
            vofref[] = nodata;
        }
    }
    restriction({c});
    boundary({vofref});//update over MPI for correct foreach
    boundary_level({vofref},max_level);
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    //Because MPI we smooth alittle differently, first we compute true vof points and normals
    //Next we will then add the smoothening to a seperate scaalar
    //Finally we can extrat the fully smooth scalar using noauto for specified region :) 
    char vofout[80];
    sprintf(vofout,"dat/vofinfo-%5.3f-p%d.dat",t,pid());
    FILE *voffile = fopen(vofout,"w");
    foreach(){
        int ref = (int)vofref[];
        struct smooths *smoothnow = vofref.smooth;
        fprintf(voffile,"%f %f %f\n",x,y,Delta);//outputs x y delta for reconstruction of vof field
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            coord n = facet_normal(point, c, s);
	        double alpha = plane_alpha(c[], n);
	        coord pc;
	        double area = plane_area_center(n, alpha, &pc);
	        if(area==0){
	            printf("Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	        }
            //set the structures point value, should apply for all
            smoothnow->points[ref][0] = (x+Delta*pc.x);
            smoothnow->points[ref][1] = (y+Delta*pc.y);
            smoothnow->points[ref][2] = n.x;
            smoothnow->points[ref][3] = n.y;
            foreach_neighbor(){
                int nref = (int)vofref[];
                if(cell.pid != pid() && c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata && smoothnow->mpicomputed[nref] < 0){
                    smoothnow->mpicomputed[nref] = cell.pid;
                    printf("%d setting request from %d @ nref:%d\n",pid(),cell.pid,nref);
                }
            }
        }
    }
    fflush(voffile);
    fclose(voffile);
    //Calculate the interface data
    //printf("(%d/%d)starting 2\n",pid(),comm_size);
    //multigrid_restriction({vofref});
    boundary({vofref});
    MPI_Barrier(MPI_COMM_WORLD);
    //After our barrier we want to gather MPI values and store them inside our structure
    //To do this we will first build a list of pid's we will need to send&recieve from, each will have to send and recieve some amount of values
    int *pidmarker = calloc(comm_size , sizeof(int));//counts wanted information
    struct smooths *smoothp = vofref.smooth;
    for(int i = 0; i < comm_size; i++){
        for(int j = 0; j < *smoothp->length; j++){
            int cur = smoothp->mpicomputed[j];
            if(cur >= 0 && cur == i){
                pidmarker[i]++;
            }
        }
        printf("%d pidmarker %d = %d\n",pid(),i,pidmarker[i]);
    }
    int **pidlocation = malloc(comm_size * sizeof(int*));//tracks location of wanted information
    for(int i = 0; i < comm_size; i++){
        pidlocation[i] = calloc(pidmarker[i] , sizeof(int));
        int pidindx = 0;
        for(int j = 0; j < *smoothp->length; j++){
            int cur = smoothp->mpicomputed[j];
            if(cur >= 0 && cur == i){
                pidlocation[i][pidindx] = j;
                printf("%d set location %d:%d = %d\n",pid(),i,pidindx,pidlocation[i][pidindx]);
                pidindx++;
            }
        }
    }
    int counts = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    for(int i = 0; i < comm_size; i++){
        //we take a 1 by 1 approach for sending and recieving, and block between each to ensure all are finished
        //if our loop is at our pid() we will send out requests for information from each PID()
        if(pid() == i){
            //If our loop is at our current pid() and not our current rank, then we recieve values
            for(int j = 0; j < comm_size; j++){
                if(pid() != j){
                    //we send the amount of expected points too all ranks except ours
                    MPI_Send(&pidmarker[j],1,MPI_INT,j,0,MPI_COMM_WORLD);
                    if(pidmarker[j] > 0){
                       //If we have points, then we also want to relay information about locations of points
                       for(int q = 0; q < pidmarker[j]; q++){
                           printf("%d sending value %d\n",pid(),pidlocation[j][q]);
                       }
                       MPI_Send(pidlocation[j],pidmarker[j],MPI_INT,j,1,MPI_COMM_WORLD);
                    }
                }
            }
            for(int j = 0; j < comm_size; j++){
                if(pid() != j && pidmarker[j] > 0){
                    //next we want to recieve our values of expected points from each rank
                    double *bufrec = malloc(nc*pidmarker[j]*sizeof(double));
                    MPI_Recv(bufrec,nc*pidmarker[j],MPI_DOUBLE,j,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    printf("%d got %d from %d\n",pid(),nc*pidmarker[j],j);
                    for(int q = 0; q < pidmarker[j]; q++){
                        smoothp->points[pidlocation[j][q]][0] = bufrec[q*nc  ];
                        smoothp->points[pidlocation[j][q]][1] = bufrec[q*nc+1];
                        smoothp->points[pidlocation[j][q]][2] = bufrec[q*nc+2];
                        smoothp->points[pidlocation[j][q]][3] = bufrec[q*nc+3];
                    }
                    free(bufrec);
                }
            }
        }
        else{
            MPI_Recv(&counts,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(counts > 0){
                printf("%d got %d from %d\n",pid(),counts,i);
                //only gather & send out values if needing too
                int *bufrec = malloc(counts * sizeof(int));
                MPI_Recv(bufrec,counts,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                for(int q = 0; q < counts; q++){
                    printf("%d got %d=%d from %d\n",pid(),q,bufrec[q],i);
                }
                //Now that we have recieved the requested locations of points, we will wrap them up and send them over
                double *bufsend = malloc(nc*counts*sizeof(double*));//we do light 'serization' really just extending out lists and rebuilding them
                for(int j = 0; j < counts; j++){
                    bufsend[j*nc  ] = smoothp->points[bufrec[j]][0];//assign x
                    bufsend[j*nc+1] = smoothp->points[bufrec[j]][1];//assign y
                    bufsend[j*nc+2] = smoothp->points[bufrec[j]][2];//assign nx
                    bufsend[j*nc+3] = smoothp->points[bufrec[j]][3];//assign ny
                }
                MPI_Send(bufsend,nc*counts,MPI_DOUBLE,i,2,MPI_COMM_WORLD);
                free(bufrec);
                free(bufsend);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);//barrier ensures all processes will stay relevant to eachother
    }
    for(int i = 0; i < comm_size; i++){
        free(pidlocation[i]);
    }
    free(pidlocation);
    free(pidmarker);
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            int ref = (int)vofref[];
            //printf("calc @ %f - level %d\n",vofref[],depth());
            struct smooths *smoothnow = vofref.smooth;
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = malloc(25*sizeof(double*));//allocated max amount of points
            for(int i = 0; i < 25; i++){
                localSpline[i] = calloc((nc/2),sizeof(double));
            }
            //First we will go though and collect all needed 
            //First in X
            int indx = 0;
            double *pnow = smoothnow->points[ref];
            foreach_neighbor(1){
                if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
                    int nref = (int)vofref[];
                    double *pnew = smoothnow->points[nref];
                    if(((pnew[2] / pnow[2]) > 0.  || (pnew[3] / pnow[3]) > 0.) && fabs(pnew[2] - pnow[2]) < 0.5 && fabs(pnew[3] - pnow[3]) < 0.5){
                        localSpline[indx][0] = pnew[0];
                        localSpline[indx][1] = pnew[1];
                        indx++;
                    }
                }
            }
            //Now we have all the points we want, so we next find our approx valuesj
            int n = 2;//order of our fit
            if(indx >= n){
                //allocate
                double *X = (double*)calloc((2*n+1) , sizeof(double));
                double *Y = (double*)calloc((n + 1) , sizeof(double));
                double **B = (double**)calloc((n+1) , sizeof(double*));
                double *A = (double*)calloc((n + 1) , sizeof(double));
                for(int i = 0; i < n+1; i++){
                    B[i] = malloc((n+2) * sizeof(double));
                }
                //calc arrays
                if(fabs(smoothnow->points[ref][3]) > fabs(smoothnow->points[ref][2])){ 
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
                    smoothnow->smoothpoints[ref][0] = smoothnow->points[ref][0];
                    smoothnow->smoothpoints[ref][1] = 0;
                    double *AP = (double*)calloc(n+1,sizeof(double));
                    for(int i = 0; i <= n; i++){
                        smoothnow->smoothpoints[ref][1] = smoothnow->smoothpoints[ref][1] + pow(smoothnow->points[ref][0],i) * A[i];
                        AP[i] = A[i] * i;
                    }
                    //and then calculate the norms using the prime
                    //First we calculate the tangent m at our point
                    double m = 0;
                    for(int i = 0; i <= n; i++){
                        if(i != 0){
                            m = m + AP[i] * pow(smoothnow->points[ref][0],i-1);
                        }
                    }
                    //normal m = -1/m
                    m = -1 * (1/m);
                    double b = (-1 * m * smoothnow->smoothpoints[ref][0]) + smoothnow->smoothpoints[ref][1];
                    //Calculate temp
                    //If were to the right side of the x center, we will calculate with x-1
                    //left side calculate with x+1?
                    double tx;
                    if(smoothnow->points[ref][2] < 0.){
                        tx = smoothnow->smoothpoints[ref][0] - 1;
                    }
                    else{
                        tx = smoothnow->smoothpoints[ref][0] + 1;
                    }
                    double ty = m * tx + b;
                    //Find direction vector to make
                    double tnormx = tx - smoothnow->smoothpoints[ref][0];
                    double tnormy = ty - smoothnow->smoothpoints[ref][1];
                    //finally we normalize
                    double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                    smoothnow->smoothpoints[ref][2] = tnormx / bottom;
                    smoothnow->smoothpoints[ref][3] = tnormy / bottom;
                    //ensure norms are correct direction before outputting
                    if(smoothnow->smoothpoints[ref][2] > 0. && !(smoothnow->points[ref][2] > 0.)){
                        smoothnow->smoothpoints[ref][2] = -1 * smoothnow->smoothpoints[ref][2];
                    }
                    else if(smoothnow->smoothpoints[ref][2] < 0. && !(smoothnow->points[ref][2] < 0.)){
                        smoothnow->smoothpoints[ref][2] = -1 * smoothnow->smoothpoints[ref][2];
                    }
                    if(smoothnow->smoothpoints[ref][3] > 0. && !(smoothnow->points[ref][3] > 0.)){
                        smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                    }
                    else if(smoothnow->smoothpoints[ref][3] < 0. && !(smoothnow->points[ref][3] < 0.)){
                        smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
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
                    smoothnow->smoothpoints[ref][1] = smoothnow->points[ref][1];
                    smoothnow->smoothpoints[ref][0] = 0;
                    double *AP = (double*)calloc(n+1,sizeof(double));
                    for(int i = 0; i <= n; i++){
                        smoothnow->smoothpoints[ref][0] = smoothnow->smoothpoints[ref][0] + pow(smoothnow->points[ref][1],i) * A[i];
                        AP[i] = A[i] * i;
                    }
                    //and then calculate the norms using the prime
                    //First we calculate the tangent m at our point
                    double m = 0;
                    for(int i = 0; i <= n; i++){
                        if(i != 0){
                            m = m + AP[i] * pow(smoothnow->points[ref][1],i-1);
                        }
                    }
                    //normal m = -1/m
                    m = -1 * (1/m);
                    double b = (-1 * m * smoothnow->smoothpoints[ref][1]) + smoothnow->smoothpoints[ref][0];
                    //Calculate temp
                    //If were to the right side of the x center, we will calculate with x-1
                    //left side calculate with x+1?
                    double ty;
                    if(smoothnow->points[ref][3] < 0.){
                        ty = smoothnow->smoothpoints[ref][1] - 1;
                    }
                    else{
                        ty = smoothnow->smoothpoints[ref][1] + 1;
                    }
                    double tx = m * ty + b;
                    //Find direction vector to make
                    double tnormx = tx - smoothnow->smoothpoints[ref][0];
                    double tnormy = ty - smoothnow->smoothpoints[ref][1];
                    //finally we normalize
                    double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                    smoothnow->smoothpoints[ref][2] = tnormx / bottom;
                    smoothnow->smoothpoints[ref][3] = tnormy / bottom;
                    //ensure normals are pointing in correct direction
                    if(smoothnow->smoothpoints[ref][2] > 0. && !(smoothnow->points[ref][2] > 0.)){
                        smoothnow->smoothpoints[ref][2] = -1 * smoothnow->smoothpoints[ref][2];
                    }
                    else if(smoothnow->smoothpoints[ref][2] < 0. && !(smoothnow->points[ref][2] < 0.)){
                        smoothnow->smoothpoints[ref][2] = -1 * smoothnow->smoothpoints[ref][2];
                    }
                    if(smoothnow->smoothpoints[ref][3] > 0. && !(smoothnow->points[ref][3] > 0.)){
                        smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                    }
                    else if(smoothnow->smoothpoints[ref][3] < 0. && !(smoothnow->points[ref][3] < 0.)){
                        smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                    }
                    free(A);
                    free(AP);
                }
                if(fabs(smoothnow->smoothpoints[ref][0] - smoothnow->points[ref][0]) > mdel/2 || fabs(smoothnow->smoothpoints[ref][0] - smoothnow->points[ref][0]) > mdel/2){
                    printf("\nStencil Error!\n");
                    printf("input point [%f,%f,%f,%f]\n",smoothnow->points[ref][0],smoothnow->points[ref][1],smoothnow->points[ref][2],smoothnow->points[ref][3]);
                    for(int i = 0; i < indx; i++){
                        printf("stencil:%d [%f,%f]\n",indx,localSpline[i][0],localSpline[i][1]);
                    }
                    printf("output point [%f,%f,%f,%f]\n",smoothnow->smoothpoints[ref][0],smoothnow->smoothpoints[ref][1],smoothnow->smoothpoints[ref][2],smoothnow->smoothpoints[ref][3]);
                }
                //freeup variables
                for(int i = 0; i < n+1; i++){
                    free(B[i]);
                }
                free(B);
                free(X);
                free(Y);
            }
            for(int i = 0;i < 25; i++){
                free(localSpline[i]);
            }
            free(localSpline);
            //if(vofref[] == 63.){
            //    double *pnt = smoothnow->smoothpoints[ref];
            //    printf("\nerrorchecker @ -%d => [%f,%f,%f,%f]\n",indx,pnt[0],pnt[1],pnt[2],pnt[3]);
            //    for(int q = 0; q < indx; q++){
            //        printf("%d - [%f,%f]\n",q,localSpline[q][0],localSpline[q][1]);
            //    }
            //    printf("\n");
            //}
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //finally we need to free the struct and its variables allocated
    free(localcalc);
}
void extract_ip_MPI(double ***parr,scalar c,scalar vofref,int *countn,int *dim,int *active_PID, double t){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    //because of weirdness with basilisk, we will write to interficial datafiles and then extract points from said files
    double ***passarr = NULL;
    scalar kappa[];//get curvature
    curvature (c, kappa);
    boundary({kappa});
    //setup MPI array
    passarr = malloc(comm_size * sizeof(double**));
    passarr[pid()] = malloc(*countn * sizeof(double*));
    int *tagarr = malloc(*countn * sizeof(int));
    for(int i = 0; i < *countn; i++){
        passarr[pid()][i] = calloc((*dim * 2) + 1,sizeof(double));
        tagarr[i] = 1;
    }
    char savename[80];
    sprintf(savename,"dat/interface-%5.3f-p%d.dat",t,pid());
    FILE *savefile = fopen(savename,"w");
    int *i = calloc(comm_size,sizeof(int));
    foreach(){
        //we go through our values and set them
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            int ref = (int)vofref[];
            struct smooths *smoothnow = vofref.smooth;
            for(int q = 0; q < *dim; q++){
                passarr[pid()][i[pid()]][q] = smoothnow->smoothpoints[ref][q];//set x
                passarr[pid()][i[pid()]][q+*dim] = smoothnow->smoothpoints[ref][q+*dim];//set xnorm
            }
            passarr[pid()][i[pid()]][*dim*2] = fabs(1./kappa[]);
            //if(isnan(passarr[pid()][i[pid()]][0]) || isnan(passarr[pid()][i[pid()]][1]) || isnan(passarr[pid()][i[pid()]][2]) || isnan(passarr[pid()][i[pid()]][3])){
            //    double *pnt = smoothnow->smoothpoints[ref];
            //    printf("error0 @ -%d => [%f,%f,%f,%f]\n",ref,pnt[0],pnt[1],pnt[2],pnt[3]);
            //} 
            i[pid()]++;
            //finally we loop through our close neighbors and grab their values
            foreach_neighbor(){
                if(c[] > 1e-6 && c[] < 1.-1e-6 && cell.pid != pid() && vofref[] != nodata){
                    int nref = (int)vofref[];
                    for(int q = 0; q < *dim; q++){
                        passarr[pid()][i[pid()]][q] = smoothnow->smoothpoints[nref][q];//set x
                        passarr[pid()][i[pid()]][q+*dim] = smoothnow->smoothpoints[nref][q+*dim];//set xnorm
                    }
                    passarr[pid()][i[pid()]][*dim*2] = fabs(1./kappa[]);
                    //if(isnan(passarr[pid()][i[pid()]][0]) || isnan(passarr[pid()][i[pid()]][1]) || isnan(passarr[pid()][i[pid()]][2]) || isnan(passarr[pid()][i[pid()]][3])){
                    //    double *pnt = smoothnow->smoothpoints[nref];
                    //    printf("error1 @ -%d => [%f,%f,%f,%f]",nref,pnt[0],pnt[1],pnt[2],pnt[3]);
                    //} 
                    i[pid()]++;
                }
            }
        }
    }
    int newcount = 0; 
    for(int q = 0; q < *countn; q++){
        //first we check our point, verify its not 0's, to ensure thge thin
        for(int p = 0; p < (*dim * 2) + 1; p++){
            if(passarr[pid()][q][p] == 0.){
                tagarr[q] = 0;
            }
            if(!tagarr[q]){
                break;//allow early break out if 2 values are 0.0; only one should ever be possibly 0. 
            }
        }
        //Next we do a backloop and verify the point hasnt been seen before
        if(tagarr[q]){
            for(int j = 0; j < q; j++){
                int counts = 0;
                for(int p = 0; p < (*dim * 2) + 1; p++){
                    //check individual values for each
                    if(passarr[pid()][q][p] == passarr[pid()][j][p]){
                        counts++;
                    }
                }
                if(counts == (*dim * 2) + 1){
                    tagarr[q] = 0;
                }
                if(!tagarr[q]){
                    break;//allow exit out
                }
            }
        }
        if(tagarr[q]){
            newcount++;
#if dimension == 2
            fprintf(savefile,"%f %f %f %f %f\n",passarr[pid()][q][0],passarr[pid()][q][1],passarr[pid()][q][2],passarr[pid()][q][3],passarr[pid()][q][4]);
#else
            fprintf(savefile,"%f %f %f %f %f %f %f\n",passarr[pid()][q][0],passarr[pid()][q][1],passarr[pid()][q][2],passarr[pid()][q][3],passarr[pid()][q][4],passarr[pid()][q][5],passarr[pid()][q][6]);
#endif
        
        }
    }
    //clean up needed
    fflush(savefile);
    fclose(savefile);
    delete({kappa});
    //Setup our array which will be ported out with the correct size :)
    double **newarr = malloc(newcount * sizeof(double*));
    for(int q = 0; q < newcount; q++){
        newarr[q] = calloc((*dim * 2) + 1,sizeof(double));
        for(int p = 0; p < (*dim * 2) + 1; p++){
            newarr[q][p] = passarr[pid()][q][p];
        }
    }
    //finally free up old arr
    for(int q = 0; q < *countn; q++){
        free(passarr[pid()][q]);
    }
    free(passarr[pid()]);
    free(passarr);
    free(tagarr);
    free(i);
    //set output vars
    *countn = newcount;
    *parr = newarr;
}
//run MPI
void calcSkeletonMPI(scalar f,double *alpha, int *dim,int max_level,double L,double t,double ***pskeleton, int *pskelelength,int *active_PID){ 
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    //printf("(%d/%d): made to start\n", pid(),comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    //for calculating local skeleton first before one full skeleton will be ran
    
    //local counts
    int *countn = calloc(comm_size,sizeof(int));
    //int thisrank = pid();
    foreach(){
        if(f[] > 1e-6 && f[] < 1-1e-6){
            //add to count from our opperating rank
            countn[pid()]++;
        }
        foreach_neighbor(){
            if(f[] > 1e-6 && f[] < 1-1e-6 && cell.pid != pid()){
                //add count from other ranks, ensures enough space, we can decrease later
                countn[pid()]++;
            }
        }
    }
    //printf("(%d/%d): counted: %d\n", pid(),comm_size,countn[pid()]);
    
    //set up for skeleton
    double mindis = 10.;
    foreach(reduction(min:mindis)){
        if(Delta < mindis){
            mindis = Delta;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("mindis:%f\n",mindis); 
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    scalar vofref[];
    vofref.refine = vofref.prolongation = smoothvof_prolongation;
    vofref.coarsen = vofref.restriction  = smoothvof_restriction;
    smooth_interface_MPI(sP,vofref,t,max_level);
    //printf("(%d/%d): smoothed\n", pid(),comm_size);
    boundary({vofref});//ensure values are updated for easy reads
    MPI_Barrier(MPI_COMM_WORLD);
    
    double **interfacePoints;
    //smooth interface points
    //next grab interface points
    extract_ip_MPI(&interfacePoints,f,vofref,&countn[pid()],dim,active_PID,t);
    delete({vofref});
    struct smooths *smooth = vofref.smooth;
    for(int i = 0; i < *smooth->length; i++){
        free(smooth->points[i]);
        free(smooth->smoothpoints[i]);
    }
    free(smooth->points);
    free(smooth->smoothpoints);
    free(smooth->mpicomputed);
    free(smooth->length);
    free(smooth);
    smooth = NULL;
    //printf("(%d/%d): got smooth\n",pid(),comm_size);
    //int holdIPM = countn[pid()];
    
    //calc skeleton
    //printf("(%d/%d): starting local skeleton\n", pid(),comm_size);
    countn[pid()] = countn[pid()] - 1;
    double skelemin = mindis * 0.1;
    char savename[80];
    double **skeleton = skeletize(interfacePoints,&countn[pid()],dim,savename,&skelemin,false,*alpha);
    
    //Next we thin out our skeleton
    //printf("(%d/%d): smoothing local skeleton\n", pid(),comm_size);
    double thindis = mindis;
    thinSkeleton(&skeleton,dim,&countn[pid()],alpha,&thindis);
    
    sprintf(savename,"dat/skeletonscatter-%5.3f-p%d.dat",t,pid());
    FILE *savefile = fopen(savename,"w");
    for(int i = 0; i < countn[pid()]; i++){
        fprintf(savefile,"%f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4]);
    } 
    fflush(savefile);
    fclose(savefile);
    

    //printf("(%d/%d): making spline\n", pid(),comm_size);
    ////create spline of skeleton
    //int skelen = 4;//n of splines
    //double del = L/pow(max_level,2);
    //double minbranchlength = 0.01;
    //int mxpt = 100; 
    //skeleReduce(skeleton,del,&minbranchlength,&countn,dim,&mxpt,t,skelen);
    //printf("(%d/%d):finished-stalling\n", pid(),comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("(%d/%d):finished-pass\n", pid(),comm_size);

    //Clean up needed
    //finally assign our needed pointers to look at the right place
    *pskelelength = countn[pid()];
    free(countn);
    *pskeleton = skeleton;
}
#endif

#endif
