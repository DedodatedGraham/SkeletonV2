/** aerobreakup of a drop 
 */
#include <time.h>
//#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
#include "two-phase-s.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "skele/skeleton.h"
#include "../colormap.h"
//#include "sandbox/lchirco/signature.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1

//draw
double max_level = 9;
double L = 12.;
double t_out = 0.001;       
double T_END = 5.0;
//double T_END = 0.5;
//double T_END = 0.001;
//double T_END = 0.1;
/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 83.33;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.01;
double uemax = 0.1;
double maxruntime = 2.5;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;


void plots(double t);

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

//signature 
scalar sign[], phii[];
int sign_lev = 7;

//skeleton tracker
scalar skelevof[];
//pid tracer
#if _MPI
int comm_size;
int *active_PID;
#endif

char restartname[80];

int main(int argc, char * argv[]){
#if _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    active_PID = calloc(comm_size , sizeof(int)); 
#endif
    //main func
    if (argc > 1){
        max_level = atof (argv[1]);
    }
    if (argc > 2){
        time_restart = atof (argv[2]);
    }
    if (argc > 3){
        T_END = atof (argv[3]);
        T_END = T_END / 1000;
        printf("T_END in => %f\n",T_END);
    }
    timerstart = clock();
    timerlast = clock();
    char dumpname[80];
    fprintf(stdout,"\n\n ERROR CHECKPOINT: %d %d \n\n",(int)time_restart,(int)(T_END * 1000));
    for(int i = (int)time_restart; i < (int)(T_END*1000); i++){
        //Here we restart and output plots
        double t = (double)i/1000.;
        sprintf(dumpname,"../dumps/dump-%5.3f",t);
        printf("%d %f\n",i,t); 
        if (!restore (file = dumpname)){
            printf("\n\n\nerror on restore @ %d\n\n\n",i);
            exit(1);
        }
        plots(t);
        run();
    }
}

/**
The initial drop is spherical. */
double x0 = 2.; double Rd = 0.5; 

char plotname[80];

void plots (double t){
    double tx = 0.;
    double vol = 0;
    foreach(reduction(+:vol) reduction(+:tx)){
        double dvf =f[]*dv();
        if(f[] > 1e-6){
            vol+= dvf;
            tx += dvf*x;
        }
    }
    tx /= vol;
    fprintf(stdout,"%f\n",tx);
    view(fov=8,tx=-(tx/L),ty=-0.5,sx=1,sy=1,camera="front",width=1000,height=1000);
    

    clear();
    
    draw_vof("f");
    squares(color = "level",min = 1, max = max_level);
    box();
    cells();
    sprintf(plotname, "Img/levels-%5.3f.png",t);
    save(file = plotname);
    
    clear();
    
    draw_vof("f");
    squares(color = "id",map = black_body,min = 0, max = 2);//0 for none, 1 for transition, 2 for skeleton
    box();
    cells();
    sprintf(plotname, "Img/svof-%5.3f.png",t);
    save(file = plotname);
    
    clear();
    
    fprintf(stdout,"done\n");
}
//TEST FUNCTIONS
//void testMPI(double t){
//    int dim = 2;
//    printf("start-%f @ %d\n",t,pid());
//    double alpha = 25 * PI / 180;/////INPUT ANGLE TO SKELETON 
//    //Calc skeleton will scan skeleton vof field and determine if new skeleton needs to be placed
//    calcSkeletonMPI(f,&alpha,&dim,max_level,L,t,active_PID);
//    printf("end-%f @ %d\n",t,pid());
//}
