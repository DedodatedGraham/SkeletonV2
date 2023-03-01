#include "skeletize.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
struct OutputXYNorm{
    scalar c;
    //FILE *fp;
    face vector s;
    int level;
};
//Method For ordering an interface, currently 2D only
double** orderInterface(double **inputPoints,struct kdleaf *kdstruct,int *length,int *dim){
    fprintf(stdout,"starting process\n");
    //Inputs are unordered points which wont work for smoothing
    double **newInterface = (double **)malloc(*length * sizeof(double*));//We attach new ordered points here
    int j = 0;
    while(j < *dim){
        newInterface[j] = (double*)malloc(sizeof(double));
        j ++;
    }
    fprintf(stdout,"alloc\n");
    double **visitedStack;//Holds the visisted points to send to kdtree
    
    //loop var
    //we will start at the first point and go on from there
    int i = 0;//index for current position of point
    double *currentpoint;
    bool runningloop = true;
    fprintf(stdout,"starting loop\n");
    while(runningloop){
        fprintf(stdout,"on point:%d \n",i);
        //running loop will complete when all points have been ordered in some way or decided to be forgotten if necessiary
        if(i == 0){
            //If first itteration then we will choose the first point, add it to stack and move on
            currentpoint = inputPoints[i];
            visitedStack[i] = inputPoints[i];
            newInterface[i][0] = inputPoints[i][0];
            newInterface[i][1] = inputPoints[i][1];
            if(*dim == 3){
                newInterface[i][2] = inputPoints[i][2];
            }
            i++;
        }
        else{
            double *lowestdistance;
            double *temppoint = getNearest(currentpoint,kdstruct,length,dim,visitedStack,&i,lowestdistance);
            if(lowestdistance != 0){
                //When we have conclusive results as there is a close point which has been uncounted
                currentpoint = temppoint;//set new value
                visitedStack[i] = temppoint;//assins old value into stack
                newInterface[i][0] = inputPoints[i][0];
                newInterface[i][1] = inputPoints[i][1];
                if(*dim == 3){
                    newInterface[i][2] = inputPoints[i][2];
                }
                i++;//increments as were moving on
            }
            else{
                fprintf(stdout,"caution empty return");
            }
        }
        if(i == *length){
            //we have counted all the points
            runningloop = false;
        }
    }
    //in future full implementation since we malloc new list to seperate dependencies, 
    //we will also free our input points as it wont get cleaned anywhere else
    //free(inputPoints)
    return newInterface;
}

void output_orderedlistxy(double **list,int *length,char path[80]){ 
    FILE * fp1 = fopen (path, "w");
    //output_facets (f,fp1);
    for(int i = 0; i < *length; i++){
        fprintf(fp1,"%g %g %d\n",list[i][0],list[i][1],i);
    }
    fflush(fp1);
    fclose(fp1);
}
void smooth(double **originalInterface,int *length, int *dim,double t){
    //For our smoothing function we are input the VOF points
    //We will first create a robust spline of our data
    //From the robust spline we will place the amount of original points all equally spaced along the spline
    //Resulting in similar quality skeletons for given data input
    //However there will be an ability to hard set the amount of points if a higher/lower resolution is desired
    //first we will create a kdtree and order the points
    
    //first create tree
    struct kdleaf *kdstruct;
    CreateStructure(originalInterface,&kdstruct,0,dim,0,*length);//make kd-struct
    fprintf(stdout,"\n");
    fprintf(stdout,"made kdtree\n");
    //order the points
    double **newpoints;
    newpoints = orderInterface(originalInterface,kdstruct,length,dim);
    fprintf(stdout,"ordered\n");
    char path[80];
    sprintf(path,"orderedinterface-%5.3f",t);
    output_orderedlistxy(newpoints,length,path);
    //Apply Spline to points
    //double **splinepoints;
    
    //clean up
    kdDestroy(kdstruct);

    //return final points
    //return splinepoints;
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

