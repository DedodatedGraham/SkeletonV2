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
#include "../colormap.h"
#include "sandbox/lchirco/signature.h"

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

void checkforThin(double t);
void runSkeleton(double t);
void plots(double t);

//signature
scalar phii[],sign[];
int sign_lev = 8;

bool hasSkeleton = false;

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
    char inname[80] = "../dumps/dump-s-10-" ;  
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
            checkforThin(ti);
            //if(hasSkeleton){
            //    runSkeleton(ti);
            //}
            plots(ti);
        }
        else{
            break;
        }
    }
}
void checkforThin(double t){
    foreach(){
        phii[] = 2*f[] - 1;
        sign[] = 7;
    }
    int l_sign = sign_lev;
    for(int ilev = depth() - 1; ilev >= l_sign; ilev--){
        foreach_level(ilev){
            if(is_refined(cell))
            restriction_average(point,phii);
        }
    }
    compute_signature_neigh_level(f,phii,sign,l_sign);
    int countnow = 0;
    //foreach_level(l_sign){
    //    //set values for each level 
    //    int val = sign[];
    //    if(val == -1 || val == 2){
    //        sign[] = -1;
    //    }
    //    else if (val == 0){
    //        sign[] = 0;
    //    }
    //    else{
    //        sign[] = 1;
    //        //here we add to our count
    //        countnow++;
    //    }
    //}
    if(countnow > 0){
        fprintf(stdout,"THIN DETECTION\n");
        hasSkeleton = true;
    }
    for(int ilev = l_sign; ilev < depth(); ilev++){
        foreach_level(ilev){
            sign.prolongation = phii.prolongation = refine_injection;
            if(is_refined(cell)){
                sign.prolongation(point,sign);
                phii.prolongation(point,phii);
            }
        }
    }
}
void runSkeleton(double t){
    fprintf(stdout,"Skeleton at %f\n",t);
    
    //setup
    char sname[80];
    sprintf (sname, "smoothskeleton-%5.3f.dat", t);
    
    double alpha = 25 * PI / 180;/////INPUT ANGLE TO SKELETON 
    double **skeleton;
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    int snr;int snd;
    double calc_time = 0;
    
    char redname[80];
    sprintf(redname, "reducedskeleton-%5.3f.dat", t);
    
    clock_t begin = clock();
    //calc minimumdistance
    double mindis = 10.;
    foreach(serial){
        if(Delta < mindis){
            mindis = Delta;
        }
    }
    fprintf(stdout,"thresh:%f\n",mindis); 
    
    
    //run smooth skeleton 
    double **sinterface = output_points_2smooth(sP,&snr,&snd,t);
    skeleton = skeletize(sinterface,&snr,&snd,sname,&mindis,false,alpha);
    
    clock_t end = clock();
    calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
    fprintf(stdout,"time took for smooth skeleton: %f\n",calc_time);
    
    //next we thin it out on alpha & mindis
    double thindis = 3*mindis;
    thinSkeleton(&skeleton,&snd,&snr,&alpha,&thindis);

    calc_time = 0;
    begin = clock();
    
    //setup spline
    int skelen = 3;//n of splines
    double del = L/pow(max_level,2);
    double minbranchlength = 0.01;
    int mxpt = 100; 
    skeleReduce(skeleton,del,&minbranchlength,&snr,&snd,&mxpt,t,skelen);
    
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

char name[80];
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
    view(fov=8,tx=-(tx/L),ty=-0.5,sx=1,sy=1,camera="front",width=2000,height=2000);
    

    clear();
    box();
    squares(color = "level",min = 1, max = max_level);
    draw_vof("f");
    sprintf(outname, "Img/levels-%5.3f.png",t);
    save(file = outname);
    
    clear();
    
    squares(color = "sign",map = black_body,linear = false,max = 2, min = -1);
    box();
    //cells(lc = {0,255,255});
    draw_vof("f");
    sprintf(name, "Img/sign-s-%d-%5.3f.png",(int)max_level,t);
    save(file = name);
    
    clear();
    
    fprintf(stdout,"done\n");
}
