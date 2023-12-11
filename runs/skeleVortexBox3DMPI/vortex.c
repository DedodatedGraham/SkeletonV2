/** aerobreakup of a drop 
 */
#include <time.h>
#include "navier-stokes/centered.h"
//#include "two-phase.h"
#include "two-phase-s.h"

#include "navier-stokes/conserving.h"
#include "tension.h"

#include "curvature.h"
#include "view.h"
#include "tag.h"
#include "skele/skeleton.h"
#include "colormap.h"
#include "curvature.h"

#define LARGE 1e36

extern scalar *interfacefid;//id 

double max_level = 6;
double t_out = 0.01;       
//double T_END = 15.;
double T_END = 2.00;

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double L = 1.;
double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 83.33;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.0001;
double uemax = 0.01;
double maxruntime = 2.5;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;

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
    sprintf(dumpname,"dumps/dump-%5.3f",0.550);
    if (!restore (file = dumpname)) {
        fraction(f,circle(x,y));
    }
}

event adapt(i++){
#if TREE
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, max_level);
#endif
}

#define velox(x,y,t,pi) ((1./pi)*cos(pi*t/50.)*sq(sin(pi*x))*2.*sin(pi*y)*cos(pi*y)*pi)
#define veloy(x,y,t,pi) ((-1./pi)*cos(pi*t/50.)*sq(sin(pi*y))*2.*sin(pi*x)*cos(pi*x)*pi)
#ifndef _PI
#define _PI 3.14159
#endif
event velocity(i++){
    //here we ensure the velocity is set
    foreach(){
        u.x[] = velox(x,y,t,_PI);
        u.y[] = veloy(x,y,t,_PI);
    }
}

char outpath[80];
event plot(t += t_out){
    view(tx=-0.5,ty=-0.5,sx=1,sy=1,camera="front",width=1000,height=1000);
    
    clear();
    sprintf(outpath,"Img/id-%5.3f.png",t);
    draw_vof("f");
    scalar fid = f.fid;
    squares(color = "fid");
    box();
    cells();
    save(file = outpath);
    clear();
    sprintf(outpath,"Img/levelint-%5.3f.png",t);
    draw_vof("f");
    squares(color = "level",min = 1, max = max_level);
    box();
    cells();
    save(file = outpath);
    
    clear();
//    
//    //sprintf(outpath,"Img/ux-%5.3f.png",t);
//    //squares(color = "u.x");
//    //box();
//    //save(file = outpath);
//    //
//    //clear();
//    //sprintf(outpath,"Img/uy-%5.3f.png",t);
//    //squares(color = "u.y");
//    //box();
//    //save(file = outpath);
//    //
//    //clear();
}

event snapshot (t += t_out; t<=T_END) {
    //Here we calculate time taken
    char dumpname[80];
    scalar * dump_list =  (scalar *){f,u,p};
#if 1
    sprintf (dumpname, "dumps/snapshot-%5.3f.gfs", t);
    output_gfs (file = dumpname, t=t,list = dump_list);
#endif

    sprintf (dumpname, "dumps/dump-%5.3f", t);
    //p.nodump = true;
    dump (file = dumpname,dump_list); // so that we can restart
    

    //finally output timer and results as this takes bulk of process
    timernow = clock();
    double cpu_time_used = ((double) (timernow - timerlast)) / (CLOCKS_PER_SEC * comm_size);
    timertotal += cpu_time_used;
    printf("calc @ t=%f\n",t);
    fprintf(stdout,"Total Time used: %fm: %fs\n",floor(timertotal / 60.),timertotal - 60 * (floor(timertotal / 60.)));
    fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
    int ncells = 0;
    foreach(reduction(+:ncells)){
        ncells++;
    }
    fprintf(stdout,"Ncells --> %d\n",ncells);
    timerlast = timernow;
}

