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
struct OutputXYZNorm{
    scalar c;
    //FILE *fp;
    face vector s;
    int level;
};
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
	            printf("Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
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
//Basilisk -> Skeleton -> Basilisk MPI comunication 
#if _MPI
void smooth_interface_MPI(struct OutputXYNorm p,vector svof, vector svofn){
    int comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    scalar c = p.c;
    vector vof[];
    vector vofn[];//temp variables used for storage
    restriction({c});
    printf("(%d,%d): we have made it\n",pid(),comm_size);
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    //Because MPI we smooth alittle differently, first we compute true vof points and normals
    //Next we will then add the smoothening to a seperate scaalar
    //Finally we can extrat the fully smooth scalar using noauto for specified region :) 
    printf("(%d/%d)starting 1\n",pid(),comm_size);
    foreach(){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            coord n = facet_normal(point, c, s);
	        double alpha = plane_alpha(c[], n);
	        coord pc;
	        double area = plane_area_center(n, alpha, &pc);
	        if(area==0){
	            printf("Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	        }
            vof.x[] = (x+Delta*pc.x);
            vof.y[] = (y+Delta*pc.y);
            vofn.x[] = n.x;
            vofn.y[] = n.y;
            //Here we capture the min and max's to prevent incorrect points from making it to the scheme when we use a grid larger than 3x3
	    }
    }
    int nc = 4;
    int grabarea = 3;//Odd number, the amount of area we want to scan around a point. Eg. if 5 we get -2,+2 from cell in each dim
    //Calculate the interface data
    printf("(%d/%d)starting 2\n",pid(),comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    boundary({vof,vofn});
    foreach(){
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
                        localSpline[indx][0] = vof.x[i,j];
                        localSpline[indx][1] = vof.y[i,j];
                        indx++;
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
            if(fabs(vofn.y[]) > fabs(vofn.x[])){ 
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
                svof.x[] = x;
                svof.y[] = 0;
                double *AP = (double*)calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    svof.y[] = svof.y[] + pow(x,i) * A[i];
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
                double b = (-1 * m * svof.x[]) + svof.y[];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double tx;
                if(vofn.x[] < 0.){
                    tx = svof.x[] - 1;
                }
                else{
                    tx = svof.x[] + 1;
                }
                double ty = m * tx + b;
                //Find direction vector to make
                double tnormx = tx - svof.x[];
                double tnormy = ty - svof.y[];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                svofn.x[] = tnormx / bottom;
                svofn.y[] = tnormy / bottom;
                //ensure norms are correct direction before outputting
                if(svofn.x[] > 0. && !(vofn.x[] > 0.)){
                    svofn.x[] = -1 * svofn.x[];
                }
                else if(svofn.x[] < 0. && !(vofn.x[] < 0.)){
                    svofn.x[] = -1 * svofn.x[];
                }
                if(svofn.y[] > 0. && !(vofn.y[] > 0.)){
                    svofn.y[] = -1 * svofn.y[];
                }
                else if(svofn.y[] < 0. && !(vofn.y[] < 0.)){
                    svofn.y[] = -1 * svofn.y[];
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
                svof.y[] = y;
                svof.x[] = 0;
                double *AP = (double*)calloc(n+1,sizeof(double));
                for(int i = 0; i <= n; i++){
                    svof.x[] = svof.x[] + pow(y,i) * A[i];
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
                double b = (-1 * m * svof.y[]) + svof.x[];
                //Calculate temp
                //If were to the right side of the x center, we will calculate with x-1
                //left side calculate with x+1?
                double ty;
                if(vofn.y[] < 0.){
                    ty = svof.y[] - 1;
                }
                else{
                    ty = svof.y[] + 1;
                }
                double tx = m * ty + b;
                //Find direction vector to make
                double tnormx = tx - svof.x[];
                double tnormy = ty - svof.y[];
                //finally we normalize
                double bottom = sqrt(pow(tnormx,2)+pow(tnormy,2));
                svofn.x[] = tnormx / bottom;
                svofn.y[] = tnormy / bottom;
                //ensure normals are pointing in correct direction
                if(svofn.x[] > 0. && !(vofn.x[] > 0.)){
                    svofn.x[] = -1 * svofn.x[];
                }
                else if(svofn.x[] < 0. && !(vofn.x[] < 0.)){
                    svofn.x[] = -1 * svofn.x[];
                }
                if(svofn.y[] > 0. && !(vofn.y[] > 0.)){
                    svofn.y[] = -1 * svofn.y[];
                }
                else if(svofn.y[] < 0. && !(vofn.y[] < 0.)){
                    svofn.y[] = -1 * svofn.y[];
                }
                free(A);
                free(AP);
            }
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
    MPI_Barrier(MPI_COMM_WORLD);
    delete({vof.x ,vof.y });//, free({vof});
    delete({vofn.x,vofn.y});//, free({vofn});
}
void extract_ip_MPI(double ***parr,scalar c,vector svof,vector svofn,int *countn,int *dim,int *active_PID, double t){
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
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            int q = 0;
            foreach_dimension(){
                passarr[pid()][i[pid()]][q] = svof.x[];//set x
                passarr[pid()][i[pid()]][q+*dim] = svofn.x[];//set xnorm
                q ++;
            }
            passarr[pid()][i[pid()]][*dim*2] = fabs(1./kappa[]);
            i[pid()]++;
            //finally we loop through our close neighbors and grab their values
            foreach_neighbor(){
                if(c[] > 1e-6 && c[] < 1.-1e-6 && cell.pid != pid()){
                    q = 0;
                    foreach_dimension(){
                        passarr[pid()][i[pid()]][q] = svof.x[];//set x
                        passarr[pid()][i[pid()]][q+*dim] = svofn.x[];//set xnorm
                        q ++;
                    }
                    passarr[pid()][i[pid()]][*dim*2] = fabs(1./kappa[]);
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
    printf("(%d/%d): made to start\n", pid(),comm_size);
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
    printf("(%d/%d): counted: %d\n", pid(),comm_size,countn[pid()]);
    
    //set up for skeleton
    double mindis = 10.;
    foreach(reduction(min:mindis)){
        if(Delta < mindis){
            mindis = Delta;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    vector svof[],svofn[];
    smooth_interface_MPI(sP,svof,svofn);
    printf("(%d/%d): smoothed\n", pid(),comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    boundary({svof,svofn});//ensure values are updated for easy reads
    
    double **interfacePoints;
    //smooth interface points
    //next grab interface points
    extract_ip_MPI(&interfacePoints,f,svof,svofn,&countn[pid()],dim,active_PID,t);
    printf("(%d/%d): got smooth\n",pid(),comm_size);
    //int holdIPM = countn[pid()];
    
    //calc skeleton
    printf("(%d/%d): starting local skeleton\n", pid(),comm_size);
    countn[pid()] = countn[pid()] - 1;
    double skelemin = mindis * 0.1;
    char savename[80];
    double **skeleton = skeletize(interfacePoints,&countn[pid()],dim,savename,&skelemin,false,*alpha);
    
    //Next we thin out our skeleton
    printf("(%d/%d): smoothing local skeleton\n", pid(),comm_size);
    double thindis = mindis * 2;
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
    printf("(%d/%d):finished-stalling\n", pid(),comm_size);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("(%d/%d):finished-pass\n", pid(),comm_size);

    //Clean up needed
    delete({svof.x ,svof.y });// free({svof });
    delete({svofn.x,svofn.y});// free({svofn});
    //finally assign our needed pointers to look at the right place
    *pskelelength = countn[pid()];
    free(countn);
    *pskeleton = skeleton;
}
#endif

#endif
