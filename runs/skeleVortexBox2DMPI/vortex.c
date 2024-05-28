/** aerobreakup of a drop 
 */
#include <time.h>
#include "navier-stokes/centered.h"
#include "two-phase.h"
//#include "two-phase-s.h"
#include "skele/skeleInclude.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
//
#include "curvature.h"
//#include "view.h"
#include "tag.h"
//#include "skele/skeleton.h"
//#include "colormap.h"
//#include "curvature.h"

#define LARGE 1e36

double max_level = 6;
double t_out = 0.01;       
//double T_END = 0.05;
double T_END = 3.0;
//double T_END = 6.00;
//double T_END = 0.11;
/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rho_L=1 [0];
double rho_V=1.297e-3 [0];
double mu_L=5.275e-3 [0];
double mu_V=9.077e-5 [0];
double femax = 0.0001;
double uemax = 0.01;
double maxruntime = 2.5;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;
double cpu_time_used = 0.;
//pid tracer
#if _MPI
int comm_size;
#endif

char restartname[80];

int main(int argc, char * argv[]){
#if _MPI
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif
    //main func
    if (argc > 1){
        max_level = atof (argv[1]);
    }
    if (argc > 2){
        time_restart = atof (argv[2]);
    }
    timerstart = clock();
    timerlast = clock();
    printf("checking %f\n",time_restart);
    origin (0,0);
    //set val
    rho1 = 1., rho2 = 1; 
    mu1 = 1, mu2 = 1;
    init_grid(1 << (int)max_level); 
    run();
}

#define circle(x,y) (sq(0.2) - (sq(x - 0.5) + sq(y - .75)))

event init (t = 0){
    char dumpname[80];
    sprintf(dumpname,"dumps/dump-%5.3f",10000.550);
    if (!restore (file = dumpname)) {
        fraction(f,circle(x,y));
    }
}

event adapt(i++){
#if TREE
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, max_level);
#endif
}
double _PI = 3.14159 [0,0];
double u0 = 1.[1,-1];
double tc = 1.[0,1];
event velocity(i++){
    //here we ensure the velocity is set
    foreach(){
        u.x[] = u0*((1./_PI)*cos(_PI*t/(tc*50.))*sq(sin(_PI*x/L0))*2.*sin(_PI*y/L0)*cos(_PI*y/L0)*_PI);
        u.y[] = u0*((-1./_PI)*cos(_PI*t/(tc*50.))*sq(sin(_PI*y/L0))*2.*sin(_PI*x/L0)*cos(_PI*x/L0)*_PI);
    }
}

event skeletizeEvent(t+=t_out){
    double **skeleton=NULL,alpha=0.1;
    int skelelength=0;
#if _MPI
    calcSkeletonMPI(f,&alpha,max_level,1.0,t,&skeleton,&skelelength,0);
#else
    calcSkeleton(f,&alpha,max_level,1.0,t,&skeleton,&skelelength,0);
#endif
    for(int q = 0; q < skelelength; q++){
        free(skeleton[q]);
    }
    free(skeleton);
}

//char outpath[80];
//event plot(t += t_out){
//    view(tx=-0.5,ty=-0.5,sx=1,sy=1,camera="front",width=1000,height=1000);
//    
//    clear();
//    sprintf(outpath,"Img/id-%5.3f.png",t);
//    draw_vof("f");
//    squares(color = "");
//    box();
//    cells();
//    save(file = outpath);
//    clear();
//    //sprintf(outpath,"Img/levelint-%5.3f.png",t);
//    //draw_vof("f");
//    //squares(color = "level",min = 1, max = max_level);
//    //box();
//    //cells();
//    //save(file = outpath);
//    //
//    //clear();
////    
////    //sprintf(outpath,"Img/ux-%5.3f.png",t);
////    //squares(color = "u.x");
////    //box();
////    //save(file = outpath);
////    //
////    //clear();
////    //sprintf(outpath,"Img/uy-%5.3f.png",t);
////    //squares(color = "u.y");
////    //box();
////    //save(file = outpath);
////    //
////    //clear();
//}

event snapshot (t += t_out; t<=T_END) {
    //Here we calculate time taken
    char dumpname[80];
    //scalar * dump_list =  (scalar *){f,u,p};

    sprintf (dumpname, "dumps/dump-%5.3f", t);
    //p.nodump = true;
    dump (file = dumpname); // so that we can restart
    

    //finally output timer and results as this takes bulk of process
    if(pid() == 0){
        timernow = clock();
        cpu_time_used = ((double) (timernow - timerlast)) / (CLOCKS_PER_SEC);
        timertotal += cpu_time_used;
        printf("calc @ t=%f - %f\n",t,max_level);
        fprintf(stdout,"Total Time used: %fm: %fs\n",floor(timertotal / 60.),timertotal - 60 * (floor(timertotal / 60.)));
        fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
        timerlast = timernow;
    }
}

