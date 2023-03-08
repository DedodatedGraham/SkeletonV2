#include "skeletize.h"
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
double* getSpline(){

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
    foreach_level_or_leaf(p.level){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            j++;
	    }
    }
    int nr = j; int nc = 4;// nc is the number of column, we initialize it with 4 because we will stor x,y norm data in those columns
    //fprintf(stdout,"Number of interfacial cells=%d\n",nr);
    int nrp = 2*nr - 1;
    *nrow = nrp;
    //fprintf(stdout,"nr=%d\n",j);
    //fprintf(stdout,"nrow=%d\n",*nrow);
    j = 0;
    double **arr = (double**)malloc((nrp+1)*sizeof(double*));
    //fprintf(stdout,"allocating %d points\n",nrp+1);
    for(int k = 0; k < nrp+1; k++){
        //fprintf(stdout,"allocating %d dim in point %d\n",nc,k);
        arr[k] = (double*)malloc(nc*sizeof(double));
    }
    //Calculate the interface data
    foreach_level_or_leaf(p.level){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            coord n = facet_normal(point, c, s);
	        double alpha = plane_alpha(c[], n);
	        coord pc;
	        double area = plane_area_center(n, alpha, &pc);
	        if(area==0){
	            fprintf(stdout,"Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	        }
	        //fprintf(stdout,"Adding in point %d\n",j);
	        arr[j][0] = x+Delta*pc.x; 
	        arr[j][1] = y+Delta*pc.y;
	        double abs = sqrt(pow(n.x,2)+pow(n.y,2));
            double tx = n.x/abs;
            double ty = n.y/abs;
            arr[j][2] = tx;  
            arr[j][3] = ty;
	        //fprintf(stdout,"Adding in point %d\n",j + 1);
	        arr[j+1][0] = (x+Delta*pc.x); 
	        arr[j+1][1] = -(y+Delta*pc.y);
            arr[j+1][2] = tx;  
            arr[j+1][3] = -ty;
            //fprintf(stdout,"loaded [%f,%f],[%f,%f]\n",arr[j][0],arr[j][1],arr[j][2],arr[j][3]);
            //fprintf(stdout,"loaded [%f,%f],[%f,%f]\n",arr[j+1][0],arr[j+1][1],arr[j+1][2],arr[j+1][3]);
	        j = j + 2;
	    }
    }
    return arr;
}

double** output_points_2smooth(struct OutputXYNorm p, int *nrow,int *ndim){
    *ndim = 2;
    scalar c = p.c;
    restriction({c});
    face vector s = p.s;
    if(!s.x.i) s.x.i = -1;
    int pj = 0;// number of interfacial cells
    scalar vofx[],vofy[];
    fprintf(stdout,"adding in:\n");
    double xmax = 0.,xmin = 0.,ymax = 0.,ymin = 0.;
    foreach_level_or_leaf(p.level){
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
            fprintf(stdout,"Point %d: [%f,%f]\n",pj,x+Delta*pc.x,y+Delta*pc.y);
            pj++;
	    }
    }
    int nr = pj; int nc = 4;// nc is the number of column, we initialize it with 4 because we will stor x,y norm data in those columns
    int nrp = 2*nr - 1;
    *nrow = nrp;
    int j = 0;
    double **arr = (double**)malloc((nrp+1)*sizeof(double*));
    for(int k = 0; k < nrp+1; k++){
        arr[k] = (double*)malloc(nc*sizeof(double));
    }
    int grabarea = 5;//Odd number, the amount of area we want to scan around a point. Eg. if 5 we get -2,+2 from cell in each dim
    //Calculate the interface data
    foreach_level_or_leaf(p.level){
        if(c[] > 1e-6 && c[] < 1.-1e-6){
            //Here we know we are currently located at an interface point we want. 
            double **localSpline = (double**)malloc((grabarea * grabarea)*sizeof(double*));//allocated max amount of points
            //First we will go though and collect all needed 
            //First in X
            int indx = 0;
            fprintf(stdout,"\ncalculating for [%f,%f]\n",x,y);
            for(int i = -1 * (grabarea - 1) / 2 ;i <= (grabarea - 1) / 2;i++){
                //Then in Y
                for(int j = - 1 * (grabarea - 1) / 2;j <= (grabarea - 1) / 2;j++){
                    //we can now go through each stencil and find all interface points along the area
                    if(c[i,j] > 1e-6 && c[i,j] < 1.-1e-6){
                        if(vofx[i,j] <= xmax && vofx[i,j] >= xmin && vofy[i,j] <= ymax && vofy[i,j] >= ymin){
                            localSpline[indx] = (double*)malloc((nc/2)*sizeof(double));
                            localSpline[indx][0] = vofx[i,j];
                            localSpline[indx][1] = vofy[i,j];
                            fprintf(stdout,"(%d,%d)captured point = [%f,%f]\n",i,j,localSpline[indx][0],localSpline[indx][1]);
                            double *splineElements = getSpline(localSpline);
                            indx++;
                            
                        }
                        else{
                            fprintf(stdout,"(%d,%d) error point = [%f,%f]\n",i,j,vofx[i,j],vofy[i,j]);
                        }
                    }
                    else{
                        fprintf(stdout,"(%d,%d) skipped point\n",i,j);
                    }
                }
            }
            //Now we have all the points we want, so we next find our spline function
            free(localSpline);
        }
    }
    return arr;
}

