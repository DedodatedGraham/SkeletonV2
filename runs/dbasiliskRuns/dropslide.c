/** aerobreakup of a drop 
 */
//#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "skele/skeleInclude.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1
double max_level = 10;
double L = 8.;
double t_out = 0.01;       
double T_END = 0.45;    


double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 78.54;        //free stream velocity
//double u0 = 107.5;
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 60;
//time var

//spline parameter
double del = 0.05/4;

double i_start = 0.;
double i_end = 39.;
//double i_start = 27.;
//double i_end = 27.;
double i_gap = 1.;
int i = 0.;
double mindis = 0.0;
int slevel = 0.;

void runSkeleton(double ti);

int main(int argc, char * argv[])
{
    //Note We scale parameters in here bc annoying comand line
    //Time Start, Time End, Lambda, Delta
    char name[80];
    char inname[80];
    if( argc > 1 ){
        sscanf(argv[1],"%lf",&i_start);
        if(i_start < 0){
            fprintf(stdout,"error corrected start");
            i_start = 0.;
        }
    }
    if( argc > 2 ){
        sscanf(argv[2],"%lf",&i_end);
        if(i_end > 39){
            fprintf(stdout,"error corrected start");
            i_end = 39.;
        }
    }
    if( argc > 3 ){
        //in del changed how small it is
        //like level so del will be 0.05 / del <- int from command line ie (1,4,100,ect.)
        sscanf(argv[3],"%lf",&del);
        del = 0.05 / del;
    }
    bool rename = false;
    if( argc > 4 ){
        sscanf(argv[4],"%s",inname);
        rename = true;
    }
    for (i = i_start; i <= i_end; i += i_gap) {
        double ti = ((double)i)/100;
        if(rename){
            char ending[10];
            sprintf(ending,"%5.3f",ti);
            char newname[80];
            strcpy(newname,inname);
            strcat(newname,ending);
            strcpy(name,newname);
            
        }
        else{
            sprintf(name, "../basiliskRuns/dump-%5.3f", ti);
        }
        fprintf (ferr, "Trying file %s, t = %f\n",name,ti);
        restore(file = name);
        fprintf (ferr, "Opening file %s, i=%d \n",name,i);
        runSkeleton(ti);
    }
}

void output_skeleinterface(char name[80],double **list,int length){
    FILE *fp = fopen(name,"w");
    for(int i = 0; i < length; i++){
        fprintf(fp,"%f %f %f %f\n",list[i][0],list[i][1],list[i][2],list[i][3]);
    }
    fflush(fp);
    fclose(fp);
}

void runSkeleton(double ti){
    //get minimumdistance
    fprintf(stdout,"\n");
    if(slevel == 0 || slevel < max_level){
        foreach(){
            if(f[] > 1e-6 && f[] < 1-1e-6){
                if(Delta < mindis || mindis == 0.){
                    mindis = Delta;
                    slevel = max_level;
                }
            }
        }
    }
    fprintf(stdout,"thresh:%f\n",mindis); 
    fprintf(stdout,"Skeleton at %f\n",ti);
    
    //setup smooth
    double alpha = 25 * PI / 180;/////INPUT ANGLE TO SKELETON 
    double **skeleton;
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    int snr;int snd;
    char sname[80];
    sprintf (sname, "smoothskeleton-%5.3f.dat", ti);
    double calc_time = 0;
    
    //run smooth
    clock_t begin = clock();
     
    double **sinterface = output_points_2smooth(sP,&snr,&snd,t);
    skeleton = skeletize(sinterface,&snr,&snd,sname,&mindis,false,alpha);
    
    clock_t end = clock();
    calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
    fprintf(stdout,"time took for smooth skeleton: %f\n",calc_time);
    
    thinSkeleton(&skeleton,&snd,&snr,&alpha);
    char redname[80];
    sprintf(redname, "reducedskeleton-%5.3f.dat", ti);
    output_skeleinterface(redname,skeleton,snr);
    calc_time = 0;
    begin = clock();
    //setup reduce
    double minbranchlength = 0.01;
    int mxpt = 100; 
    //run reduce
    //double del = 0.05;
    //double del = 0.05/4;
    int n = 3;
    skeleReduce(skeleton,del,&minbranchlength,&snr,&snd,&mxpt,ti,n);
    
    end = clock();
    calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
    fprintf(stdout,"time took for reduced skeleton: %f\n",calc_time);
    //cleanup data
    for(int i = 0; i < snr;i++){
        free(skeleton[i]);
    }
    free(skeleton);
    skeleton = NULL;

    fflush(ferr);
    fprintf(stdout,"\n");
}



