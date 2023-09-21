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
#include "../colormap.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1

//draw
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

//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

//dump var
double i_start = 0.;
double i_end = 2000.;
double time_restart = 0.;
double i_gap = 1;
int i = 0.;

void plots(double t);

int main(int argc, char * argv[]){
    //main func

    if (argc > 1){
        max_level = atof (argv[1]);
    }
    if (argc > 2){
        time_restart = atof (argv[2]);
    }
    if(time_restart != 0.){
        i_start = time_restart;
    }
    timerstart = clock();
    timerlast = clock();
    char name[80];
    char inname[80] = "../dumps/dump-ns-10-" ;  
    for (i = i_start; i <= i_end; i += i_gap) {
        double ti = ((double)i)/1000;
        char ending[10];
        sprintf(ending,"%5.3f",ti);
        char newname[80];
        strcpy(newname,inname);
        strcat(newname,ending);
        strcpy(name,newname);
        fprintf (ferr, "Trying file %s, t = %f\n",name,ti);
        bool go = restore(file = name);
        if(go){
            fprintf (ferr, "Opening file %s, i=%d \n",name,i);
            plots(ti);
        }
        else{
            break;
        }
    }
}
char outname[80];
void plots (double t){
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
    sprintf(outname, "Img/levels-%5.3f.png",t);
    save(file = outname);
    clear();
    fprintf(stdout,"done\n");
}
