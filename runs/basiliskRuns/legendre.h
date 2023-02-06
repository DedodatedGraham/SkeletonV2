#include <stdio.h>
#include <math.h>
#include <time.h>

#define PI 3.1416

float LegPol(int n, float x){
/*   if(x< -1 || x>1){
        return 0;
    }
    else{*/
    if(n == 0){
	return 1.;
    }
    else if(n == 1){
        return x;
    }
    else if(n == 2){
	return (1./2.)*(3*pow(x,2) - 1);
    }
    else if(n == 3){
	return (1.0/2.0)*(5*pow(x,3) - 3*x);
    }
    else if(n == 4){
        return (1./8.)*(35*pow(x,4) - 30*pow(x,2) + 3);
    }
    else if(n == 5){
	return (1./8.)*(63*pow(x,5) - 70*pow(x,3) + 15*x);
    }
    else if(n == 6){
	return (1./16.)*(231*pow(x,6) - 315*pow(x,4) + 105*pow(x,2) - 5);
    }
    else if(n == 7){
        return (1./16.)*(429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) - 35*x);
    }
    else if(n == 8){
        return (1./128.)*(6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260*pow(x,2) + 35);
    }
    else if(n == 9){
        return (1./128.)*(12155*pow(x,9) - 25740*pow(x,7) + 18018*pow(x,5) - 4620*pow(x,3) + 315*x);
    }
    else if(n == 10){
        return (1./256.)*(46189*pow(x,10) - 109396*pow(x,8) + 90090*pow(x,6) - 30030*pow(x,4) + 3465*pow(x,2) - 63);
    }
 // }
return 0;
}


void Swap(double *a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int Partition(double **arr, int low, int high, int index){
    double pivot = arr[high][index]; //initializing 1st pivot, index is the index of the column whose elements are used as sorting key
    int i = low - 1; // intializing second pointer

    for(int j = low; j < high; j++){
        if(arr[j][index] <= pivot){
	    i++;
	    Swap(&arr[i][index], &arr[j][index]);// swap the current element with the second pointer
	    Swap(&arr[i][index-1], &arr[j][index-1]);// swap the current element with the second pointer
	    Swap(&arr[i][index-2], &arr[j][index-2]);// swap the current element with the second pointer
	}
    }
    Swap(&arr[i+1][index], &arr[high][index]);// swap the pivot with the second pointer at the end of whole pass
    Swap(&arr[i+1][index-1], &arr[high][index-1]);// swap the pivot with the second pointer at the end of whole pass
    Swap(&arr[i+1][index-2], &arr[high][index-2]);// swap the pivot with the second pointer at the end of whole pass
    int i_pivot = i+1; //this is the new pivot location
    //display(arr);
return i_pivot;
}

void QuickSort(double **arr, int low, int high, int index){
    if(low < high){
        double piv = Partition(arr,low,high,index);
//        printf("QuickSort\n");   
        QuickSort(arr, low, piv-1, index);
        QuickSort(arr, piv+1, high, index);
//    printf("low = %d high = %d\n",low, high);
    }
}


//Extract the interfacial points. here we are extracting the center of the cut surface
struct OutputXYTheta{
    scalar c;
    //FILE *fp;
    face vector s;
    int level;
};

double** output_points_xytheta(struct OutputXYTheta p, double xc, int *nrow){
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
    int nr = j; int nc = 3;// nc is the number of column, we initialize it with 3 because we will stor x,y theta data in those columns
    //fprintf(stdout,"Number of interfacial cells=%d\n",nr);
    *nrow = j;
    //fprintf(p.fp,"#x y\n");
    j = 0;
    //Lets allocate memory to store the interface data as an array
    //fprintf(stdout,"DEBUG: Lets allocate memory\n");
    double **arr = (double**)malloc(nr*sizeof(double*));
    for(int k = 0; k < nr; k++){
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
	    //fprintf(p.fp, "%g %g\n",x+Delta*pc.x, y+Delta*pc.y);
	    arr[j][0] = x+Delta*pc.x; 
	    arr[j][1] = y+Delta*pc.y;
            arr[j][2] = (arr[j][0] - xc < 0) ? (180 + atan( arr[j][1]/(arr[j][0]-xc+0.000000001))*180/PI) : (atan(arr[j][1]/(arr[j][0]-xc+0.00000001))*180/PI);
	    j++;
	}
    }
    //fprintf(stdout,"DEBUG:Number of interfacial cells=%d\n",j);
    QuickSort(arr, 0, nr-1, nc-1);
    //fprintf(stdout,"DEBUG: Size of arr=%g\n",sizeof(**arr));
    //Safely break program if extracted data is not enough to get correct Mode Coefficient C.
    if(arr[0][2] > 1. || arr[nr-1][2] < 179.){
	//fprintf(stdout,"ALERT!!: Missing out some interfacial points. Try adjusting mesh level\n");
        //fprintf(stdout,"theta_min = %g, theta_max = %g\n",arr[0][2], arr[nr-1][2]);
	/*//Increase maxlevel if some points are missing
	p.c = c; p.s = s; p.level = p.level+1;
	int Level = p.level;
	fprintf(stdout,"maxlevel = %d\n",p.level);
	fprintf(stdout,"xc = %g\n",xc);
	output_points_xytheta(p,xc,&nr);
	if(Level>16){
	    exit(0);
	}*/
    }

return arr;
    
free(arr);
    //arr = NULL;
    //fflush(p.fp);
}

// Print an intedded array
void Display(double **array, int r, int c){
   for(int i = 0; i < r; i++){
       for(int j = 0; j < c; j++){
           fprintf(stdout,"%0.3f ", array[i][j]);
       }
       printf("\n");
   }
}
//Calculate Fourier-Legendre coefficient 
double** fLegCoeff(double **arr, int n_max, int nr, double xc, double R_avg){
    //Allocate memory and initialize modes, say n=0 to 10
    int * n_arr = (int*)malloc((n_max+1)*sizeof(int));
    for(int i = 0; i <= n_max; i++){
        n_arr[i] = i;
	//fprintf(stdout,"DEBUG: n=%d\n",n_arr[i]);
    }
    
    //fprintf(stdout,"xc=%g\n",xc);
    
    //Calculate R_avg from interface data. However for even mode R_avg calculation from interface is little faulty. So we are doing it from vof
    //double R = 0.; 
    double *r = (double*)malloc(nr*sizeof(double));
    for(int j = 0; j < nr; j++){
        r[j] = sqrt( pow(arr[j][0]-xc , 2) + pow(arr[j][1] , 2) );
	//fprintf(stdout,"r = %g\n",r[j]);
	//R = R+r[j];
    }
    //double R_avg = R/nr;
   
    //fprintf(stdout,"R_avg=%g\n",R_avg);
    
    //Allocate memory for storing mode coeff C[i]
    double ** C = (double**)malloc((n_max+1)*sizeof(double*));
    for(int i = 0 ; i <= n_max; i++){
	C[i] = (double*)malloc(2*sizeof(double));// We are going to store mode number n and coresponding coeff c in the array C. so there are 2 cols
    }

    double dtheta;
    for(int i = 0; i <= n_max; i++){
	double c = 0;
           	
	// Integrate using trapezoidal rule
	for(int j = 1; j < nr; j++){
	    dtheta = fabs( cos(arr[j][2]*PI/180) - cos(arr[j-1][2]*PI/180) );
	    c = c + 0.5*( (r[j]/R_avg - 1)*LegPol(n_arr[i], cos(arr[j][2]*PI/180)) + (r[j-1]/R_avg -1)*LegPol(n_arr[i], cos(arr[j-1][2]*PI/180)) ) * dtheta;
	}
	
	/*
	// Integrate with naive sum method
        for(int j = 0; j < nr; j++){
	    if(j==0){
	        dtheta = 0.;
	        //dtheta = fabs( cos(arr[j+1][2]*PI/180) - cos(arr[j][2]*PI/180) ); 
	    }
	    else{
	        dtheta = fabs( cos(arr[j][2]*PI/180) - cos(arr[j-1][2]*PI/180) );
	    }
	    //Calculate Coeff
	    //printf("DEBUG:i=%d n=%d\n",i,n_arr[i]);
	    //c = c + (sqrt( pow(arr[j][0]-xc, 2) + pow(arr[j][1], 2) )/R_avg - 1) * LegPol(n_arr[i], cos(arr[j][2]*PI/180)) * dtheta;
	    c = c + (r[j]/R_avg - 1) * LegPol(n_arr[i], cos(arr[j][2]*PI/180)) * dtheta;
	}*/
             	
        C[i][0] = n_arr[i];
        C[i][1] = c*(2*n_arr[i]+1)/2;
        //printf("n=%0.0f,C=%g\n",C[i][0],C[i][1]);
        //fprintf(fp,"%d %0.3f\n",i,C[i]);
    }
return C;
free(n_arr);
free(r);
free(C);
}

