#include "skeletize.h"
//Extract the interfacial points. here we are extracting the center of the cut surface
struct OutputXYNorm{
    scalar c;
    //FILE *fp;
    face vector s;
    int level;
};

double** output_points_xynorm(struct OutputXYTheta p, int *nrow,int *ndim){
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
	        arr[j+1][0] = -(x+Delta*pc.x); 
	        arr[j+1][1] = y+Delta*pc.y;
            arr[j+1][2] = -tx;  
            arr[j+1][3] = ty;
            //fprintf(stdout,"loaded [%f,%f],[%f,%f]\n",arr[j][0],arr[j][1],arr[j][2],arr[j][3]);
            //fprintf(stdout,"loaded [%f,%f],[%f,%f]\n",arr[j+1][0],arr[j+1][1],arr[j+1][2],arr[j+1][3]);
	        j = j + 2;
	    }
    }
    return arr;
}
