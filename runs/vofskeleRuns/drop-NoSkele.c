/** aerobreakup of a drop 
 */
#include <time.h>
//#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "skele/skeleton.h"
//#include "sandbox/lchirco/signature.h"
#include "colormap.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1

//draw
scalar plotS[];
scalar s[],phii[];
double max_level = 10;
double L = 12.;
double t_out = 0.001;       
double T_END = 1.50;    
//double T_END = 0.1;
/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 87.8;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 2.5;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);


int main(int argc, char * argv[]){
    //main func
    if (argc > 1){
        max_level = atof (argv[1]);
    }
    timerstart = clock();
    timerlast = clock();
    char name[80];
    if(time_restart == 0.0){
        size (L);
        origin(0,0);
        init_grid (128);
    }
    else{
        sprintf(name, "dump-%5.3f", time_restart);
        fprintf (ferr, "Trying file %s, t = %f\n",name,time_restart);
        restore(file = name);
    }

    /**
    The dimensionless parameters can then be computed based on 
    rho_l, sigma, and D. */
    rho1 = 1., rho2 = rho_V; 
    mu1 = mu_L, mu2 = mu_V;
    f.sigma = 1.;
 
    
    //TOLERANCE = 1.e-5; 
    //NITERMIN = 5;
    run();
}

/**
The initial drop is spherical. */
double x0 = 2.; double Rd = 0.5; 


event init (t = 0){
    if (!restore (file = "dump")) {
        refine (  sq(y-(L/2))+sq(x-x0) < sq(1.2*Rd) && sq(y-(L/2))+sq(x-x0) > sq(0.8*Rd) && level < max_level);
        fraction (f, -sq(y-(L/2))-sq(x-x0)+sq(Rd));

        /** It is important to initialize the flow field. If u.x=0, Poisson
        solver crahses. If u.x=u0*(1-f[]), the interface velocity is too large
        and will induce erroneous deformation at early time and wrong pressure 
        field. As along as we use u.x=a*u0*(1-f[]), with a<0.01, then the results 
        look normal and not sensitive to a. */
        foreach() {
            u.x[]=1.e-7*u0*(1.-f[]);
        }
        boundary ({f,u.x});
    }
}
event snapshot (t += t_out; t<=T_END) {
    //Here we calculate time taken
    char name[80];
#if 1
    sprintf (name, "dumps/snapshot-ns-%d-%5.3f.gfs",(int)max_level, t);
    output_gfs (file = name, t=t, list = {f,u,p});
#endif

    sprintf (name, "dumps/dump-ns-%d-%5.3f",(int)max_level, t);
    p.nodump = false;
    dump (file = name); // so that we can restart

    //finally output timer and results as this takes bulk of process
    timernow = clock();
    double cpu_time_used = ((double) (timernow - timerlast)) / (CLOCKS_PER_SEC * 28);
    timertotal += cpu_time_used;
    fprintf(stdout,"Total Time used: %fm: %fs\n",floor(timertotal / 60.),timertotal - 60 * (floor(timertotal / 60.)));
    fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
    timerlast = timernow;
}
char name[80];
event plots (t += t_out; t<=T_END){
    fprintf(stdout,"plotting\n");
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
    view(fov=10,tx=-(tx/L),ty=-0.5,sx=1,sy=1,camera="front",width=1000,height=1000);
    

    clear();
    box();
    squares(color = "level",min = 1, max = max_level);
    draw_vof("f");
    sprintf(name, "Img/levels-ns-%d-%5.3f.png",(int)max_level,t);
    save(file = name);
    clear();
    fprintf(stdout,"done\n");
}

/**
Adapt mesh based on the volume fraction. 
*/

event adapt (i=10;i++) {
#if 0
  vector q[];
  foreach()
    foreach_dimension()
      q.x[] = u.x[]*(rho1*f[]+rho2*(1.-f[]));
  boundary ((scalar *){q});
  adapt_wavelet ({f,q}, (double[]){femax,uemax,uemax}, minlevel = 3, maxlevel = max_level);
#else
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, maxlevel = max_level);
  //adapt_wavelet2((scalar *){f,u.x,u.y},(double []){femax,uemax,uemax},(int []){max_level, max_level-2, max_level-2},3);
#endif 
}

