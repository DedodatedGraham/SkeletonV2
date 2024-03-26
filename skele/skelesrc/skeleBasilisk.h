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
    if(getDistance(p1,p2) > neededdis){
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
            if(ppoint[2] < maxz && ppoint[2] > minz)return 1;
#else
            return 1;
#endif
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
static void unsmoothvof_prolongation(Point point, scalar s){
    double val = s[];
    if(s[] != nodata){
        double *ppoint = s.smooth->smoothpoints[(int)s[]];
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
int checknorm(double *p0,double *p1){
    double dp = 0.;
    for(int i = 0; i < dimension; i++){
        dp = dp + p0[dimension+i]  * p1[dimension+i];
    }
    if(dp >= 0){
        return 1;
    }
    return 0;
}

#if _MPI
void unsmooth_interface_MPI(struct OutputXYNorm p,scalar vofref,double t,int max_level){
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
        if(i <= pid() - 1){
            lstart += localcalc[i];
        }
        *calcl += localcalc[i];
    }
    smooth->length = calcl;
    double **smoothpoints = malloc(*calcl * sizeof(double*));
    int *beencalc = calloc(*calcl , sizeof(int));
    for(int i = 0; i < *calcl; i++){
        smoothpoints[i] = calloc(nc,sizeof(double));
        beencalc[i] = -1; 
    }
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
    //restriction({c});
    //boundary({vofref});//update over MPI for correct foreach
    //restriction({vofref});
    //for(int q = max_level - 1; q >= 0; q--){
    //    boundary_level({vofref},q);
    //}
    //boundary({vofref});//update over MPI for correct foreach
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    //Because MPI we smooth alittle differently, first we compute true vof points and normals
    //Next we will then add the smoothening to a seperate scaalar
    //Finally we can extrat the fully smooth scalar using noauto for specified region :) 
    foreach(){
        int ref = (int)vofref[];
        struct smooths *smoothnow = vofref.smooth;
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            coord n = facet_normal(point, c, s);
	        double aalpha = plane_alpha(c[], n);
            normalize(&n);
#if dimension == 2
	        coord pc;
	        plane_area_center(n, aalpha, &pc);
            smoothnow->smoothpoints[ref][0] = (x+Delta*pc.x);
            smoothnow->smoothpoints[ref][1] = (y+Delta*pc.y);
            smoothnow->smoothpoints[ref][2] = n.x;
            smoothnow->smoothpoints[ref][3] = n.y;
#else
	        coord pc[12];
	        int m = facets(n, aalpha, pc, 1.);
            for(int i = 0; i < m; i++){
                smoothnow->smoothpoints[ref][0] = smoothnow->smoothpoints[ref][0] + (x+Delta*pc[i].x);
                smoothnow->smoothpoints[ref][1] = smoothnow->smoothpoints[ref][1] + (y+Delta*pc[i].y);
                smoothnow->smoothpoints[ref][2] = smoothnow->smoothpoints[ref][2] + (z+Delta*pc[i].z);
            }
            smoothnow->smoothpoints[ref][0] = smoothnow->smoothpoints[ref][0] / m;
            smoothnow->smoothpoints[ref][1] = smoothnow->smoothpoints[ref][1] / m;
            smoothnow->smoothpoints[ref][2] = smoothnow->smoothpoints[ref][2] / m;
            smoothnow->smoothpoints[ref][3] = n.x;
            smoothnow->smoothpoints[ref][4] = n.y;
            smoothnow->smoothpoints[ref][5] = n.z;
#endif
        }
    }
    //Calculate the interface data
    //multigrid_restriction({vofref});
    //boundary({vofref});
    MPI_Barrier(MPI_COMM_WORLD);
    foreach(){
        struct smooths *smoothnow = vofref.smooth;
        foreach_neighbor(){
            int nref = (int)vofref[];
            if(cell.pid != pid() && c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata && smoothnow->mpicomputed[nref] < 0 && !smoothnow->smoothpoints[nref][0]){
                smoothnow->mpicomputed[nref] = cell.pid;
            }
        }
    }
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
    }
    int **pidlocation = malloc(comm_size * sizeof(int*));//tracks location of wanted information
    for(int i = 0; i < comm_size; i++){
        pidlocation[i] = calloc(pidmarker[i] , sizeof(int));
        int pidindx = 0;
        for(int j = 0; j < *smoothp->length; j++){
            int cur = smoothp->mpicomputed[j];
            if(cur >= 0 && cur == i){
                pidlocation[i][pidindx] = j;
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
                       MPI_Send(pidlocation[j],pidmarker[j],MPI_INT,j,1,MPI_COMM_WORLD);
                    }
                }
            }
            for(int j = 0; j < comm_size; j++){
                if(pid() != j && pidmarker[j] > 0){
                    //next we want to recieve our values of expected points from each rank
                    double *bufrec = malloc(nc*pidmarker[j]*sizeof(double));
                    MPI_Recv(bufrec,nc*pidmarker[j],MPI_DOUBLE,j,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    for(int q = 0; q < pidmarker[j]; q++){
                        smoothp->smoothpoints[pidlocation[j][q]][0] = bufrec[q*nc  ];
                        smoothp->smoothpoints[pidlocation[j][q]][1] = bufrec[q*nc+1];
                        smoothp->smoothpoints[pidlocation[j][q]][2] = bufrec[q*nc+2];
                        smoothp->smoothpoints[pidlocation[j][q]][3] = bufrec[q*nc+3];
#if dimension == 3
                        smoothp->smoothpoints[pidlocation[j][q]][4] = bufrec[q*nc+4];
                        smoothp->smoothpoints[pidlocation[j][q]][5] = bufrec[q*nc+5];
#endif
                    }
                    free(bufrec);
                }
            }
        }
        else{
            MPI_Recv(&counts,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(counts > 0){
                //only gather & send out values if needing too
                int *bufrec = malloc(counts * sizeof(int));
                MPI_Recv(bufrec,counts,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                //Now that we have recieved the requested locations of points, we will wrap them up and send them over
                double *bufsend = malloc(nc*counts*sizeof(double*));//we do light 'serization' really just extending out lists and rebuilding them
                for(int j = 0; j < counts; j++){
                    bufsend[j*nc  ] = smoothp->smoothpoints[bufrec[j]][0];
                    bufsend[j*nc+1] = smoothp->smoothpoints[bufrec[j]][1];
                    bufsend[j*nc+2] = smoothp->smoothpoints[bufrec[j]][2];
                    bufsend[j*nc+3] = smoothp->smoothpoints[bufrec[j]][3];
#if dimension == 3
                    bufsend[j*nc+4] = smoothp->smoothpoints[bufrec[j]][4];
                    bufsend[j*nc+5] = smoothp->smoothpoints[bufrec[j]][5];
#endif
                }
                MPI_Send(bufsend,nc*counts,MPI_DOUBLE,i,2,MPI_COMM_WORLD);
                free(bufrec);
                free(bufsend);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);//barrier ensures all processes will stay relevant to eachother
    }
    vofref.smooth = smoothp;
    for(int i = 0; i < comm_size; i++){
        free(pidlocation[i]);
    }
    free(pidlocation);
    free(pidmarker);
    free(localcalc);
    MPI_Barrier(MPI_COMM_WORLD);
}
void smooth_interface_MPI(struct OutputXYNorm p,scalar vofref,double t,int max_level){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    int nc = dimension*2;
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
        if(i <= pid() - 1){
            lstart += localcalc[i];
        }
        *calcl += localcalc[i];
    }
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
    //restriction({vofref});
    for(int q = max_level - 1; q >= 0; q--){
        boundary_level({vofref},q);
    }
    boundary({vofref});//update over MPI for correct foreach
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    //Because MPI we smooth alittle differently, first we compute true vof points and normals
    //Next we will then add the smoothening to a seperate scaalar
    //Finally we can extrat the fully smooth scalar using noauto for specified region :) 
    foreach(){
        int ref = (int)vofref[];
        struct smooths *smoothnow = vofref.smooth;
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            coord n = facet_normal(point, c, s);
	        double aalpha = plane_alpha(c[], n);
            normalize(&n);
#if dimension == 2
	        coord pc;
	        plane_area_center(n, aalpha, &pc);
            smoothnow->points[ref][0] = (x+Delta*pc.x);
            smoothnow->points[ref][1] = (y+Delta*pc.y);
            smoothnow->points[ref][2] = n.x;
            smoothnow->points[ref][3] = n.y;
#else
	        coord pc[12];
	        int m = facets(n, aalpha, pc, 1.);
            for(int i = 0; i < m; i++){
                smoothnow->points[ref][0] = smoothnow->points[ref][0] + (x+Delta*pc[i].x);
                smoothnow->points[ref][1] = smoothnow->points[ref][1] + (y+Delta*pc[i].y);
                smoothnow->points[ref][2] = smoothnow->points[ref][2] + (z+Delta*pc[i].z);
            }
            smoothnow->points[ref][0] = smoothnow->points[ref][0] / m;
            smoothnow->points[ref][1] = smoothnow->points[ref][1] / m;
            smoothnow->points[ref][2] = smoothnow->points[ref][2] / m;
            smoothnow->points[ref][3] = n.x;
            smoothnow->points[ref][4] = n.y;
            smoothnow->points[ref][5] = n.z;
#endif
        }
    }
    //Calculate the interface data
    //multigrid_restriction({vofref});
    boundary({vofref});
    MPI_Barrier(MPI_COMM_WORLD);
    foreach(){
        struct smooths *smoothnow = vofref.smooth;
        foreach_neighbor(){
            int nref = (int)vofref[];
            if(cell.pid != pid() && c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata && smoothnow->mpicomputed[nref] < 0 && !smoothnow->points[nref][0]){
                smoothnow->mpicomputed[nref] = cell.pid;
            }
        }
    }
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
    }
    int **pidlocation = malloc(comm_size * sizeof(int*));//tracks location of wanted information
    for(int i = 0; i < comm_size; i++){
        pidlocation[i] = calloc(pidmarker[i] , sizeof(int));
        int pidindx = 0;
        for(int j = 0; j < *smoothp->length; j++){
            int cur = smoothp->mpicomputed[j];
            if(cur >= 0 && cur == i){
                pidlocation[i][pidindx] = j;
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
                       MPI_Send(pidlocation[j],pidmarker[j],MPI_INT,j,1,MPI_COMM_WORLD);
                    }
                }
            }
            for(int j = 0; j < comm_size; j++){
                if(pid() != j && pidmarker[j] > 0){
                    //next we want to recieve our values of expected points from each rank
                    double *bufrec = malloc(nc*pidmarker[j]*sizeof(double));
                    MPI_Recv(bufrec,nc*pidmarker[j],MPI_DOUBLE,j,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    for(int q = 0; q < pidmarker[j]; q++){
                        smoothp->points[pidlocation[j][q]][0] = bufrec[q*nc  ];
                        smoothp->points[pidlocation[j][q]][1] = bufrec[q*nc+1];
                        smoothp->points[pidlocation[j][q]][2] = bufrec[q*nc+2];
                        smoothp->points[pidlocation[j][q]][3] = bufrec[q*nc+3];
#if dimension == 3
                        smoothp->points[pidlocation[j][q]][4] = bufrec[q*nc+4];
                        smoothp->points[pidlocation[j][q]][5] = bufrec[q*nc+5];
#endif
                    }
                    free(bufrec);
                }
            }
        }
        else{
            MPI_Recv(&counts,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(counts > 0){
                //only gather & send out values if needing too
                int *bufrec = malloc(counts * sizeof(int));
                MPI_Recv(bufrec,counts,MPI_INT,i,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                //Now that we have recieved the requested locations of points, we will wrap them up and send them over
                double *bufsend = malloc(nc*counts*sizeof(double*));//we do light 'serization' really just extending out lists and rebuilding them
                for(int j = 0; j < counts; j++){
                    bufsend[j*nc  ] = smoothp->points[bufrec[j]][0];
                    bufsend[j*nc+1] = smoothp->points[bufrec[j]][1];
                    bufsend[j*nc+2] = smoothp->points[bufrec[j]][2];
                    bufsend[j*nc+3] = smoothp->points[bufrec[j]][3];
#if dimension == 3
                    bufsend[j*nc+4] = smoothp->points[bufrec[j]][4];
                    bufsend[j*nc+5] = smoothp->points[bufrec[j]][5];
#endif
                }
                MPI_Send(bufsend,nc*counts,MPI_DOUBLE,i,2,MPI_COMM_WORLD);
                free(bufrec);
                free(bufsend);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);//barrier ensures all processes will stay relevant to eachother
    }
    vofref.smooth = smoothp;
    for(int i = 0; i < comm_size; i++){
        free(pidlocation[i]);
    }
    free(pidlocation);
    free(pidmarker);
    struct smooths *smoothnow = vofref.smooth;
    int nls = pow(5,dimension);
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            int ref = (int)vofref[];
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = malloc(nls*sizeof(double*));//allocated max amount of points
            for(int i = 0; i < nls; i++){
                localSpline[i] = calloc((nc/2),sizeof(double));
            }
            //First we will go though and collect all needed 
            //First in X
            int indx = 0;
            double *pnow = smoothnow->points[ref];
            foreach_neighbor(){
                if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
                    int nref = (int)vofref[];
                    //printf("nref->%d/%f\n",nref,vofref[]);
                    double *pnew = smoothnow->points[nref];
                    if(pnew[0] != 0. && pnew[1] != 0. && checknorm(pnow,pnew)){
                        localSpline[indx][0] = pnew[0];
                        localSpline[indx][1] = pnew[1];
#if dimension == 3
                        localSpline[indx][2] = pnew[2];
#endif
                        indx++;
                    }
                }
            }
            //Now we have all the points we want, so we next find our approx valuesj
            int n = 2;//order of our fit
#if dimension == 2
            //2D smoothing
            if(indx > n){
                //allocate
                double *X = calloc((2*n+1) , sizeof(double));
                double *Y = calloc((n + 1) , sizeof(double));
                double **B= malloc((n+1)   * sizeof(double*));
                double *A = calloc((n + 1) , sizeof(double));
                for(int i = 0; i < n+1; i++){
                    B[i] = malloc((n+2) * sizeof(double));
                }
                //calc arrays
                int jx,jy;
                if(fabs(smoothnow->points[ref][3]) > fabs(smoothnow->points[ref][2])){ 
                    jx=0,jy=1;
                }
                else{
                    jx=1,jy=0;
                        
                }
                for(int i = 0; i <= 2*n; i++){
                    X[i] = 0;
                    for(int j = 0; j < indx; j++){
                        X[i] = X[i] + pow(localSpline[j][jx],i);
                    }
                }
                for(int i = 0; i <= n; i++){
                    Y[i] = 0;
                    for(int j = 0; j < indx; j ++){
                        Y[i] = Y[i] + pow(localSpline[j][jx],i)*localSpline[j][jy];
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
                smoothnow->smoothpoints[ref][jx] = smoothnow->points[ref][jx];
                smoothnow->smoothpoints[ref][jy] = 0;
                double *AP = calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    smoothnow->smoothpoints[ref][jy] = smoothnow->smoothpoints[ref][jy] + pow(smoothnow->points[ref][jx],i) * A[i];
                    AP[i] = A[i] * i;
                }
                //and then calculate the norms using the prime
                //First we calculate the tangent m at our point
                double m = 0;
                for(int i = 0; i <= n; i++){
                    if(i != 0){
                        m = m + AP[i] * pow(smoothnow->points[ref][jx],i-1);
                    }
                }
                //normal m = -1/m
                m = -1 * (1/m);
                double b = (-1 * m * smoothnow->smoothpoints[ref][jx]) + smoothnow->smoothpoints[ref][jy];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double tx;
                if(smoothnow->points[ref][jx+dimension] < 0.){
                    tx = smoothnow->smoothpoints[ref][jx] - 1;
                }
                else{
                    tx = smoothnow->smoothpoints[ref][jx] + 1;
                }
                double ty = m * tx + b;
                //Find direction vector to make
                double tnormx = tx - smoothnow->smoothpoints[ref][jx];
                double tnormy = ty - smoothnow->smoothpoints[ref][jy];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                smoothnow->smoothpoints[ref][jx+dimension] = tnormx / bottom;
                smoothnow->smoothpoints[ref][jy+dimension] = tnormy / bottom;
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
                //freeup variables
                for(int i = 0; i < n+1; i++){
                    free(B[i]);
                }
                free(B);
                free(X);
                free(Y);
            }
#else
            n = 2;//in 3D we settle for a plane normally
            //3D smoothening using z(x,y) = a0 + a1x + .. + anx^n + an+1y + ... + a2ny^n
            if(indx > n){
                //allocate
                double **B = malloc((2*n+1) * sizeof(double*));
                double  *A = calloc((2*n+1) , sizeof(double ));
                for(int i = 0; i < 2*n+1; i++){
                    B[i]   = calloc((2*n+2) , sizeof(double ));
                }
                //calc arrays
                int jx,jy,jz;
                if(fabs(smoothnow->points[ref][5]) >= fabs(smoothnow->points[ref][3]) && fabs(smoothnow->points[ref][5]) >= fabs(smoothnow->points[ref][4])){ 
                    jx=0,jy=1,jz=2;
                }
                else if(fabs(smoothnow->points[ref][4]) >= fabs(smoothnow->points[ref][3]) && fabs(smoothnow->points[ref][4]) >= fabs(smoothnow->points[ref][5])){ 
                    jx=0,jy=2,jz=1;
                }
                else{
                    jx=1,jy=2,jz=0;
                }
                //Set X, half x half y
                //Finally set B
                //X*X terms
                for(int i = 0; i < 2*n + 1; i++){
                    //shift between x & y values
                    int indi=jx,powi=i;
                    if(i > n)indi=jy,powi=i-n;
                    for(int j = 0; j < 2*n + 1; j++){
                        int indj=jx,powj=j;
                        if(j > n)indj=jy,powj=j-n;
                        for(int q = 0; q < indx; q++){
                            B[i][j] = B[i][j] + pow(localSpline[q][indi],powi) * pow(localSpline[q][indj],powj);
                        }
                    }
                    for(int q = 0; q < indx; q++){
                        B[i][2*n+1] = B[i][2*n+1] + localSpline[q][jz] * pow(localSpline[q][indi],powi);
                    }
                }
                getCoeffGE(2*n+1,2*(n+1),&B,&A);
                //Finally we will Get our current point
                smoothnow->smoothpoints[ref][jx] = smoothnow->points[ref][jx];
                smoothnow->smoothpoints[ref][jy] = smoothnow->points[ref][jy];
                smoothnow->smoothpoints[ref][jz] = 0;
                for(int i = 0; i < 2*n + 1; i++){
                    int addi = jx,powi = i; 
                    if(i > n)addi = jy,powi = i-n;
                    double addup = A[i] * pow(smoothnow->smoothpoints[ref][addi],powi);
                    smoothnow->smoothpoints[ref][jz] = smoothnow->smoothpoints[ref][jz] + addup;
                }
                double dzdx = 0.,dzdy = 0.;
                for(int i = 1; i <= n; i++){
                    dzdx = dzdx + A[i  ] * i * pow(smoothnow->smoothpoints[ref][jx],i-1);
                    dzdy = dzdy + A[i+n] * i * pow(smoothnow->smoothpoints[ref][jy],i-1);
                }
                double tnormx = dzdx;
                double tnormy = dzdy;
                double tnormz = -1;
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2)+pow(tnormz,2));
                smoothnow->smoothpoints[ref][jx+dimension] = tnormx / bottom;
                smoothnow->smoothpoints[ref][jy+dimension] = tnormy / bottom;
                smoothnow->smoothpoints[ref][jz+dimension] = tnormz / bottom;
                //ensure norms are correct direction before outputting
                if(smoothnow->smoothpoints[ref][3] > 0. && !(smoothnow->points[ref][3] > 0.)){
                    smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                }
                else if(smoothnow->smoothpoints[ref][3] < 0. && !(smoothnow->points[ref][3] < 0.)){
                    smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                }
                if(smoothnow->smoothpoints[ref][4] > 0. && !(smoothnow->points[ref][4] > 0.)){
                    smoothnow->smoothpoints[ref][4] = -1 * smoothnow->smoothpoints[ref][4];
                }
                else if(smoothnow->smoothpoints[ref][4] < 0. && !(smoothnow->points[ref][4] < 0.)){
                    smoothnow->smoothpoints[ref][4] = -1 * smoothnow->smoothpoints[ref][4];
                }
                if(smoothnow->smoothpoints[ref][5] > 0. && !(smoothnow->points[ref][5] > 0.)){
                    smoothnow->smoothpoints[ref][5] = -1 * smoothnow->smoothpoints[ref][5];
                }
                else if(smoothnow->smoothpoints[ref][5] < 0. && !(smoothnow->points[ref][5] < 0.)){
                    smoothnow->smoothpoints[ref][5] = -1 * smoothnow->smoothpoints[ref][5];
                }
                free(A);
                //freeup variables
                for(int i = 0; i < 2*n + 1; i++){
                    free(B[i]);
                }
                free(B);
            }
#endif
            else{
                for(int q = 0; q < dimension * 2; q++)smoothnow->smoothpoints[ref][q] = smoothnow->points[ref][q];
            }
            for(int i = 0;i < nls; i++){
                free(localSpline[i]);
            }
            free(localSpline);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    vofref.smooth = smooth;//sets the field to point at our structure
    //finally we need to free the struct and its variables allocated
    free(localcalc);
}
void extract_ip_MPI(double ***parr,scalar c,scalar vofref,int *countn, double t){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    //because of weirdness with basilisk, we will write to interficial datafiles and then extract points from said files
    double ***passarr = NULL;
    scalar kappa[];//get curvature
    curvature (c, kappa);
    boundary({kappa});
    int lcount = 0;
    foreach(noauto){
        //Count local counts
        if(c[] > 1e-6 &&  c[] < 1-1e-6){
            lcount++;
        }
    }
    *countn = lcount;
    //setup MPI array
    passarr = malloc(comm_size * sizeof(double**));
    passarr[pid()] = malloc(*countn * sizeof(double*));
    int *tagarr = malloc(*countn * sizeof(int));
    for(int i = 0; i < *countn; i++){
        passarr[pid()][i] = calloc((dimension * 2) + 1,sizeof(double));
        tagarr[i] = 1;
    }
    char savename[80];
    sprintf(savename,"dat/interface-%5.3f-p%03d.txt",t,pid());
    FILE *savefile = fopen(savename,"w");
    int *i = calloc(comm_size,sizeof(int));
    foreach(){
        //we go through our values and set them
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            int ref = (int)vofref[];
            struct smooths *smoothnow = vofref.smooth;
            for(int q = 0; q < dimension; q++){
                passarr[pid()][i[pid()]][q] = smoothnow->smoothpoints[ref][q];//set x
                passarr[pid()][i[pid()]][q+dimension] = smoothnow->smoothpoints[ref][q+dimension];//set xnorm
            }
            passarr[pid()][i[pid()]][dimension*2] = fabs(1./kappa[]);
            i[pid()]++;
        }
    }
    int newcount = 0; 
    for(int q = 0; q < *countn; q++){
        //first we check our point, verify its not 0's, to ensure thge thin
        for(int p = 0; p < (dimension * 2) + 1; p++){
            if(passarr[pid()][q][p] == 0.){
                tagarr[q] = 0;
            }
            if(!tagarr[q]){
                break;//allow early break out if 2 values are 0.0; only one should ever be possibly 0. 
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
    int q = 0;
    for(int t = 0; t < *countn; t++){
        if(tagarr[t]){// && passarr[pid()][t][*dim] != 0. && passarr[pid()][t][*dim+1] != 0.){
            newarr[q] = calloc((dimension * 2) + 1,sizeof(double));
            for(int p = 0; p < (dimension * 2) + 1; p++){
                newarr[q][p] = passarr[pid()][t][p];
            }
            //printf("newinterface Skeleton - %d => [%f,%f,%f,%f]\n",q,passarr[pid()][q][0],passarr[pid()][q][1],passarr[pid()][q][2],passarr[pid()][q][3]);
            q++;
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
void calcSkeletonMPI(scalar f,double *alpha,int max_level,double L,double t,double ***pskeleton, int *pskelelength,int doSmooth){ 
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    //set up for skeleton
    int countn = 0;
    double mindis = 10.;
    foreach(reduction(min:mindis)){
        if(Delta < mindis){
            mindis = Delta;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    face vector ps = {{-1}};
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;sP.s = ps;
    scalar vofref[];
    if(doSmooth){
        vofref.refine = vofref.prolongation = smoothvof_prolongation;
    }
    else{
        vofref.refine = vofref.prolongation = unsmoothvof_prolongation;
    }
    vofref.coarsen = vofref.restriction  = smoothvof_restriction;
    foreach(){
        vofref[] = -1.;
    }
    if(doSmooth){
        smooth_interface_MPI(sP,vofref,t,max_level);
        boundary({vofref});//ensure values are updated for easy reads
    }
    else{
        unsmooth_interface_MPI(sP,vofref,t,max_level);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    double **interfacePoints;
    //smooth interface points
    //next grab interface points
    extract_ip_MPI(&interfacePoints,f,vofref,&countn,t);
    delete({vofref});
    struct smooths *smooth = vofref.smooth;
    for(int i = 0; i < *smooth->length; i++){
        if(doSmooth)free(smooth->points[i]);
        free(smooth->smoothpoints[i]);
    }
    if(doSmooth)free(smooth->points);
    free(smooth->smoothpoints);
    free(smooth->mpicomputed);
    free(smooth->length);
    free(smooth);
    smooth = NULL;
    double skelemin = mindis * 0.5;
    char savename[80];
    double **skeleton = NULL; 
    printf("counted: %d\n",countn);
    skeletizeMPI(interfacePoints,&countn,savename,&skelemin,*alpha,&skeleton);
    MPI_Barrier(MPI_COMM_WORLD);
    //Next we thin out our skeleton
    double thindis = mindis;//when r is more than delta/2 :(
    char noname[80];
    sprintf(noname,"dat/skeletonthin-%5.3f-p%03d.txt",t,pid());
    thinSkeleton(&skeleton,&countn,alpha,&thindis,noname);
    MPI_Barrier(MPI_COMM_WORLD);
    
    sprintf(savename,"dat/skeletonscatter-%5.3f-p%03d.txt",t,pid());
    FILE *savefile = fopen(savename,"w");
    for(int i = 0; i < countn; i++){
#if dimension == 2
        fprintf(savefile,"%f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4]);
#else
        fprintf(savefile,"%f %f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4],skeleton[i][5]);
#endif
    } 
    fflush(savefile);
    fclose(savefile);
    

    ////create spline of skeleton
    //int skelen = 4;//n of splines
    //double del = L/pow(max_level,2);
    //double minbranchlength = 0.01;
    //int mxpt = 100; 
    //skeleReduce(skeleton,del,&minbranchlength,&countn,dim,&mxpt,t,skelen);
    MPI_Barrier(MPI_COMM_WORLD);

    //Clean up needed
    //finally assign our needed pointers to look at the right place
    *pskelelength = countn;
    *pskeleton = skeleton;
}
#else
//serial calculation
void unsmooth_interface(struct OutputXYNorm p,scalar vofref,double t,int max_level){
#if dimension == 2
    int nc = 4;
#else
    int nc = 6;
#endif
    scalar c = p.c;
    struct smooths *smooth = malloc(sizeof(struct smooths));
    //next we calc our total id across all MPI, to ensure we build our ID's correctly
    int localcalc = 0;
    //this will make life vastly easier when smoothening across boundaries
    double mdel = HUGE;
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            localcalc++;
            if(Delta < mdel){
                mdel = Delta;
            }
        }
    }
    int *sendl = malloc(sizeof(int));
    *sendl = localcalc;
    smooth->length = sendl;
    double **smoothpoints = malloc(localcalc * sizeof(double*));
    for(int i = 0; i < localcalc; i++){
        smoothpoints[i] = calloc(nc,sizeof(double));
    }
    smooth->smoothpoints = smoothpoints;
    vofref.smooth = smooth;//sets the field to point at our structure
    //so now we have allocated a structure belonging to our ref which has enough size to hold all of the points, regardless of if they cross a boundary
    //however we will only assign in points which fall across our foreach_neighbor; first we allocate our ID's
    int lstart = 0;
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            vofref[] = lstart++; 
        }
        else{
            vofref[] = nodata;
        }
    }
    //restriction({c});
    //boundary({vofref});
    //for(int q = max_level - 1; q >= 0; q--){
    //    boundary_level({vofref},q);
    //}
    //boundary({vofref});
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    foreach(){
        int ref = (int)vofref[];
        struct smooths *smoothnow = vofref.smooth;
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            coord n = facet_normal(point, c, s);
	        double aalpha = plane_alpha(c[], n);
            normalize(&n);
#if dimension == 2
	        coord pc;
	        plane_area_center(n, aalpha, &pc);
            smoothnow->smoothpoints[ref][0] = (x+Delta*pc.x);
            smoothnow->smoothpoints[ref][1] = (y+Delta*pc.y);
            smoothnow->smoothpoints[ref][2] = n.x;
            smoothnow->smoothpoints[ref][3] = n.y;
#else
	        coord pc[12];
	        int m = facets(n, aalpha, pc, 1.);
            for(int i = 0; i < m; i++){
                smoothnow->smoothpoints[ref][0] = smoothnow->points[ref][0] + (x+Delta*pc[i].x);
                smoothnow->smoothpoints[ref][1] = smoothnow->points[ref][1] + (y+Delta*pc[i].y);
                smoothnow->smoothpoints[ref][2] = smoothnow->points[ref][2] + (z+Delta*pc[i].z);
            }
            smoothnow->smoothpoints[ref][0] = smoothnow->points[ref][0] / m;
            smoothnow->smoothpoints[ref][1] = smoothnow->points[ref][1] / m;
            smoothnow->smoothpoints[ref][2] = smoothnow->points[ref][2] / m;
            smoothnow->smoothpoints[ref][3] = n.x;
            smoothnow->smoothpoints[ref][4] = n.y;
            smoothnow->smoothpoints[ref][5] = n.z;
#endif
        }
    }
    //boundary({vofref});
}

void smooth_interface(struct OutputXYNorm p,scalar vofref,double t,int max_level){
#if dimension == 2
    int nc = 4;
#else
    int nc = 6;
#endif
    scalar c = p.c;
    struct smooths *smooth = malloc(sizeof(struct smooths));
    //next we calc our total id across all MPI, to ensure we build our ID's correctly
    int localcalc = 0;
    //this will make life vastly easier when smoothening across boundaries
    double mdel = HUGE;
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            localcalc++;
            if(Delta < mdel){
                mdel = Delta;
            }
        }
    }
    int *sendl = malloc(sizeof(int));
    *sendl = localcalc;
    smooth->length = sendl;
    double **vofpoints = malloc(localcalc * sizeof(double*));
    double **smoothpoints = malloc(localcalc * sizeof(double*));
    //int *beencalc = calloc(localcalc , sizeof(int));
    for(int i = 0; i < localcalc; i++){
        vofpoints[i] = calloc(nc,sizeof(double));
        smoothpoints[i] = calloc(nc,sizeof(double));
        //beencalc[i] = -1; 
    }
    smooth->points = vofpoints;
    smooth->smoothpoints = smoothpoints;
    //smooth->mpicomputed = beencalc;
    vofref.smooth = smooth;//sets the field to point at our structure
    //so now we have allocated a structure belonging to our ref which has enough size to hold all of the points, regardless of if they cross a boundary
    //however we will only assign in points which fall across our foreach_neighbor; first we allocate our ID's
    int lstart = 0;
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            vofref[] = lstart++; 
        }
        else{
            vofref[] = nodata;
        }
    }
    restriction({c});
    boundary({vofref});
    for(int q = max_level - 1; q >= 0; q--){
        boundary_level({vofref},q);
    }
    boundary({vofref});
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    foreach(){
        int ref = (int)vofref[];
        struct smooths *smoothnow = vofref.smooth;
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            coord n = facet_normal(point, c, s);
	        double aalpha = plane_alpha(c[], n);
            normalize(&n);
#if dimension == 2
	        coord pc;
	        plane_area_center(n, aalpha, &pc);
            smoothnow->points[ref][0] = (x+Delta*pc.x);
            smoothnow->points[ref][1] = (y+Delta*pc.y);
            smoothnow->points[ref][2] = n.x;
            smoothnow->points[ref][3] = n.y;
#else
	        coord pc[12];
	        int m = facets(n, aalpha, pc, 1.);
            for(int i = 0; i < m; i++){
                smoothnow->points[ref][0] = smoothnow->points[ref][0] + (x+Delta*pc[i].x);
                smoothnow->points[ref][1] = smoothnow->points[ref][1] + (y+Delta*pc[i].y);
                smoothnow->points[ref][2] = smoothnow->points[ref][2] + (z+Delta*pc[i].z);
            }
            smoothnow->points[ref][0] = smoothnow->points[ref][0] / m;
            smoothnow->points[ref][1] = smoothnow->points[ref][1] / m;
            smoothnow->points[ref][2] = smoothnow->points[ref][2] / m;
            smoothnow->points[ref][3] = n.x;
            smoothnow->points[ref][4] = n.y;
            smoothnow->points[ref][5] = n.z;
#endif
        }
    }
    //Calculate the interface data
    boundary({vofref});
    struct smooths *smoothnow = vofref.smooth;
    int nls = pow(5,dimension);
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            int ref = (int)vofref[];
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = malloc(nls*sizeof(double*));//allocated max amount of points
            for(int i = 0; i < nls; i++){
                localSpline[i] = calloc((nc/2),sizeof(double));
            }
            //First we will go though and collect all needed 
            //First in X
            int indx = 0;
            double *pnow = smoothnow->points[ref];
            foreach_neighbor(){
                if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
                    int nref = (int)vofref[];
                    double *pnew = smoothnow->points[nref];
                    if(pnew[0] != 0. && pnew[1] != 0. && checknorm(pnow,pnew)){
                        localSpline[indx][0] = pnew[0];
                        localSpline[indx][1] = pnew[1];
#if dimension == 3
                        localSpline[indx][2] = pnew[2];
#endif
                        indx++;
                    }
                }
            }
            //Now we have all the points we want, so we next find our approx valuesj
            int n = 2;//order of our fit
#if dimension == 2
            //2D smoothing
            if(indx > n){
                //allocate
                double *X = calloc((2*n+1) , sizeof(double));
                double *Y = calloc((n + 1) , sizeof(double));
                double **B= malloc((n+1)   * sizeof(double*));
                double *A = calloc((n + 1) , sizeof(double));
                for(int i = 0; i < n+1; i++){
                    B[i] = malloc((n+2) * sizeof(double));
                }
                //calc arrays
                int jx,jy;
                if(fabs(smoothnow->points[ref][3]) > fabs(smoothnow->points[ref][2])){ 
                    jx=0,jy=1;
                }
                else{
                    jx=1,jy=0;
                        
                }
                for(int i = 0; i <= 2*n; i++){
                    X[i] = 0;
                    for(int j = 0; j < indx; j++){
                        X[i] = X[i] + pow(localSpline[j][jx],i);
                    }
                }
                for(int i = 0; i <= n; i++){
                    Y[i] = 0;
                    for(int j = 0; j < indx; j ++){
                        Y[i] = Y[i] + pow(localSpline[j][jx],i)*localSpline[j][jy];
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
                smoothnow->smoothpoints[ref][jx] = smoothnow->points[ref][jx];
                smoothnow->smoothpoints[ref][jy] = 0;
                double *AP = calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    smoothnow->smoothpoints[ref][jy] = smoothnow->smoothpoints[ref][jy] + pow(smoothnow->points[ref][jx],i) * A[i];
                    AP[i] = A[i] * i;
                }
                //and then calculate the norms using the prime
                //First we calculate the tangent m at our point
                double m = 0;
                for(int i = 0; i <= n; i++){
                    if(i != 0){
                        m = m + AP[i] * pow(smoothnow->points[ref][jx],i-1);
                    }
                }
                //normal m = -1/m
                m = -1 * (1/m);
                double b = (-1 * m * smoothnow->smoothpoints[ref][jx]) + smoothnow->smoothpoints[ref][jy];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double tx;
                if(smoothnow->points[ref][jx+dimension] < 0.){
                    tx = smoothnow->smoothpoints[ref][jx] - 1;
                }
                else{
                    tx = smoothnow->smoothpoints[ref][jx] + 1;
                }
                double ty = m * tx + b;
                //Find direction vector to make
                double tnormx = tx - smoothnow->smoothpoints[ref][jx];
                double tnormy = ty - smoothnow->smoothpoints[ref][jy];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                smoothnow->smoothpoints[ref][jx+dimension] = tnormx / bottom;
                smoothnow->smoothpoints[ref][jy+dimension] = tnormy / bottom;
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
                //freeup variables
                for(int i = 0; i < n+1; i++){
                    free(B[i]);
                }
                free(B);
                free(X);
                free(Y);
            }
#else
            n = 2;//in 3D we settle for a plane normally
            //3D smoothening using z(x,y) = a0 + a1x + .. + anx^n + an+1y + ... + a2ny^n
            if(indx > n){
                //allocate
                double **B = malloc((2*n+1) * sizeof(double*));
                double  *A = calloc((2*n+1) , sizeof(double ));
                for(int i = 0; i < 2*n+1; i++){
                    B[i]   = calloc((2*n+2) , sizeof(double ));
                }
                //calc arrays
                int jx,jy,jz;
                if(fabs(smoothnow->points[ref][5]) >= fabs(smoothnow->points[ref][3]) && fabs(smoothnow->points[ref][5]) >= fabs(smoothnow->points[ref][4])){ 
                    jx=0,jy=1,jz=2;
                }
                else if(fabs(smoothnow->points[ref][4]) >= fabs(smoothnow->points[ref][3]) && fabs(smoothnow->points[ref][4]) >= fabs(smoothnow->points[ref][5])){ 
                    jx=0,jy=2,jz=1;
                }
                else{
                    jx=1,jy=2,jz=0;
                }
                //Set X, half x half y
                //Finally set B
                //X*X terms
                for(int i = 0; i < 2*n + 1; i++){
                    //shift between x & y values
                    int indi=jx,powi=i;
                    if(i > n)indi=jy,powi=i-n;
                    for(int j = 0; j < 2*n + 1; j++){
                        int indj=jx,powj=j;
                        if(j > n)indj=jy,powj=j-n;
                        for(int q = 0; q < indx; q++){
                            B[i][j] = B[i][j] + pow(localSpline[q][indi],powi) * pow(localSpline[q][indj],powj);
                        }
                    }
                    for(int q = 0; q < indx; q++){
                        B[i][2*n+1] = B[i][2*n+1] + localSpline[q][jz] * pow(localSpline[q][indi],powi);
                    }
                }
                getCoeffGE(2*n+1,2*(n+1),&B,&A);
                //Finally we will Get our current point
                smoothnow->smoothpoints[ref][jx] = smoothnow->points[ref][jx];
                smoothnow->smoothpoints[ref][jy] = smoothnow->points[ref][jy];
                smoothnow->smoothpoints[ref][jz] = 0;
                for(int i = 0; i < 2*n + 1; i++){
                    int addi = jx,powi = i; 
                    if(i > n)addi = jy,powi = i-n;
                    double addup = A[i] * pow(smoothnow->smoothpoints[ref][addi],powi);
                    smoothnow->smoothpoints[ref][jz] = smoothnow->smoothpoints[ref][jz] + addup;
                }
                double dzdx = 0.,dzdy = 0.;
                for(int i = 1; i <= n; i++){
                    dzdx = dzdx + A[i  ] * i * pow(smoothnow->smoothpoints[ref][jx],i-1);
                    dzdy = dzdy + A[i+n] * i * pow(smoothnow->smoothpoints[ref][jy],i-1);
                }
                double tnormx = dzdx;
                double tnormy = dzdy;
                double tnormz = -1;
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2)+pow(tnormz,2));
                smoothnow->smoothpoints[ref][jx+dimension] = tnormx / bottom;
                smoothnow->smoothpoints[ref][jy+dimension] = tnormy / bottom;
                smoothnow->smoothpoints[ref][jz+dimension] = tnormz / bottom;
                //ensure norms are correct direction before outputting
                if(smoothnow->smoothpoints[ref][3] > 0. && !(smoothnow->points[ref][3] > 0.)){
                    smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                }
                else if(smoothnow->smoothpoints[ref][3] < 0. && !(smoothnow->points[ref][3] < 0.)){
                    smoothnow->smoothpoints[ref][3] = -1 * smoothnow->smoothpoints[ref][3];
                }
                if(smoothnow->smoothpoints[ref][4] > 0. && !(smoothnow->points[ref][4] > 0.)){
                    smoothnow->smoothpoints[ref][4] = -1 * smoothnow->smoothpoints[ref][4];
                }
                else if(smoothnow->smoothpoints[ref][4] < 0. && !(smoothnow->points[ref][4] < 0.)){
                    smoothnow->smoothpoints[ref][4] = -1 * smoothnow->smoothpoints[ref][4];
                }
                if(smoothnow->smoothpoints[ref][5] > 0. && !(smoothnow->points[ref][5] > 0.)){
                    smoothnow->smoothpoints[ref][5] = -1 * smoothnow->smoothpoints[ref][5];
                }
                else if(smoothnow->smoothpoints[ref][5] < 0. && !(smoothnow->points[ref][5] < 0.)){
                    smoothnow->smoothpoints[ref][5] = -1 * smoothnow->smoothpoints[ref][5];
                }
                free(A);
                //freeup variables
                for(int i = 0; i < 2*n + 1; i++){
                    free(B[i]);
                }
                free(B);
            }
#endif
            else{
                for(int q = 0; q < dimension * 2; q++)smoothnow->smoothpoints[ref][q] = smoothnow->points[ref][q];
            }
            for(int i = 0;i < nls; i++){
                free(localSpline[i]);
            }
            free(localSpline);
        }
    }
    vofref.smooth = smooth;//sets the field to point at our structure
    //finally we need to free the struct and its variables allocated
}
void extract_ip(double ***parr,scalar c,scalar vofref,int *countn, double t){
    //because of weirdness with basilisk, we will write to interficial datafiles and then extract points from said files
    double **passarr = NULL;
    scalar kappa[];//get curvature
    curvature (c, kappa);
    boundary({kappa});
    foreach(){
        //Count local counts
        if(c[] > 1e-6 &&  c[] < 1-1e-6){
            *countn = *countn+1;
        }
    }
    //setup MPI array
    passarr = malloc(*countn * sizeof(double*));
    int *tagarr = malloc(*countn * sizeof(int));
    for(int i = 0; i < *countn; i++){
        passarr[i] = calloc((dimension * 2) + 1,sizeof(double));
        tagarr[i] = 1;
    }
    char savename[80];
    sprintf(savename,"dat/interface-%5.3f-p%03d.txt",t,0);
    FILE *savefile = fopen(savename,"w");
    int i = 0;
    foreach(){
        //we go through our values and set them
        if(c[] > 1e-6 && c[] < 1.-1e-6 && vofref[] != nodata){
            //printf("(%d)\n",i);
            int ref = (int)vofref[];
            struct smooths *smoothnow = vofref.smooth;
            for(int q = 0; q < dimension; q++){
                passarr[i][q] = smoothnow->smoothpoints[ref][q];//set x
                passarr[i][q+dimension] = smoothnow->smoothpoints[ref][q+dimension];//set xnorm
            }
            passarr[i][dimension*2] = fabs(1./kappa[]);
            i++;
            //finally we loop through our close neighbors and grab their values
        }
    }
    int newcount = 0; 
    for(int q = 0; q < *countn; q++){
        //printf("(%d/%d)\n",q,*countn);
        //first we check our point, verify its not 0's, to ensure thge thin
        for(int p = 0; p < (dimension * 2) + 1; p++){
            if(passarr[q][p] == 0.){
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
                for(int p = 0; p < (dimension * 2) + 1; p++){
                    //check individual values for each
                    if(passarr[q][p] == passarr[j][p]){
                        counts++;
                    }
                }
                if(counts == (dimension * 2) + 1){
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
            fprintf(savefile,"%f %f %f %f %f\n",passarr[q][0],passarr[q][1],passarr[q][2],passarr[q][3],passarr[q][4]);
#else
            fprintf(savefile,"%f %f %f %f %f %f %f\n",passarr[q][0],passarr[q][1],passarr[q][2],passarr[q][3],passarr[q][4],passarr[q][5],passarr[q][6]);
#endif   
        }
    }
    //clean up needed
    fflush(savefile);
    fclose(savefile);
    delete({kappa});
    //Setup our array which will be ported out with the correct size :)
    double **newarr = malloc(newcount * sizeof(double*));
    int q = 0;
    for(int t = 0; t < *countn; t++){
        if(tagarr[t]){
            newarr[q] = calloc((dimension * 2) + 1,sizeof(double));
            for(int p = 0; p < (dimension * 2) + 1; p++){
                newarr[q][p] = passarr[t][p];
            }
            q++;
        }
    }
    //finally free up old arr
    for(int q = 0; q < *countn; q++){
        free(passarr[q]);
    }
    free(passarr);
    free(tagarr);
    //set output vars
    *countn = newcount;
    *parr = newarr;
}
void calcSkeleton(scalar f,double *alpha,int max_level,double L,double t,double ***pskeleton, int *pskelelength,int doSmooth){ 
    //set up for skeleton
    int countn = 0;
    double mindis = 10.;
    foreach(){
        if(Delta < mindis){
            mindis = Delta;
        }
    }
    face vector ps = {{-1}};
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;sP.s = ps;
    scalar vofref[];
    if(doSmooth){
        vofref.refine = vofref.prolongation = smoothvof_prolongation;
    }
    else{
        vofref.refine = vofref.prolongation = unsmoothvof_prolongation;
    }
    vofref.coarsen = vofref.restriction  = smoothvof_restriction;
    foreach(){
        vofref[] = -1.;
    }
    //printf("setup\n");
    if(doSmooth){
        smooth_interface(sP,vofref,t,max_level);
        boundary({vofref});//ensure values are updated for easy reads
    }
    else{
        unsmooth_interface(sP,vofref,t,max_level);
    }
    //printf("smooth\n");
    
    double **interfacePoints;
    //smooth interface points
    //next grab interface points
    extract_ip(&interfacePoints,f,vofref,&countn,t);
    //printf("extract\n");
    delete({vofref});
    struct smooths *smooth = vofref.smooth;
    for(int i = 0; i < *smooth->length; i++){
        if(doSmooth)free(smooth->points[i]);
        free(smooth->smoothpoints[i]);
    }
    if(doSmooth)free(smooth->points);
    free(smooth->smoothpoints);
    free(smooth->length);
    free(smooth);
    smooth = NULL;
    
    //calc skeleton
    double skelemin = mindis;
    char savename[80];
    double **skeleton = NULL; 
    skeletize(interfacePoints,&countn,savename,&skelemin,*alpha,&skeleton);
    //printf("skele\n");
    //Next we thin out our skeleton
    double thindis = mindis;//when r is more than delta/2 :(
    char noname[80];
    sprintf(noname,"dat/skeletonthin-%5.3f.txt",t);
    thinSkeleton(&skeleton,&countn,alpha,&thindis,noname);
    sprintf(savename,"dat/skeletonscatter-%5.3f.txt",t);
    FILE *savefile = fopen(savename,"w");
    for(int i = 0; i < countn; i++){
#if dimension == 2
        fprintf(savefile,"%f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4]);
#else
        fprintf(savefile,"%f %f %f %f %f %f\n",skeleton[i][0],skeleton[i][1],skeleton[i][2],skeleton[i][3],skeleton[i][4],skeleton[i][5]);
#endif
    } 
    //printf("save\n");
    fflush(savefile);
    fclose(savefile);
    

    ////create spline of skeleton
    //int skelen = 4;//n of splines
    //double del = L/pow(max_level,2);
    //double minbranchlength = 0.01;
    //int mxpt = 100; 
    //skeleReduce(skeleton,del,&minbranchlength,&countn,dim,&mxpt,t,skelen);

    //Clean up needed
    //finally assign our needed pointers to look at the right place
    *pskelelength = countn;
    *pskeleton = skeleton;
}
#endif

#endif
