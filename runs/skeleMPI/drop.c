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
#include "colormap.h"
#include "sandbox/lchirco/signature.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1

//draw
double max_level = 10;
double L = 12.;
double t_out = 0.001;       
double T_END = 0.005;    
//double T_END = 0.1;
/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 78.54;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.01;
double maxruntime = 2.5;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;
double timertotal = 0.;

bool hasSkeleton = false;

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

//signature 
scalar sign[], phii[];
int sign_lev = 8;
int main(int argc, char * argv[]){
    //main func
    if (argc > 1){
        max_level = atof (argv[1]);
    }
    if (argc > 2){
        time_restart = atof (argv[2]);
    }
    timerstart = clock();
    timerlast = clock();
    char name[80];
    if(time_restart == 0.0){
        size (L);
        origin(0,0);
        init_grid (128);

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
    else{
        time_restart /= 1000;
        sprintf(name, "dump-s-%2.0f-%5.3f", max_level,time_restart);
        fprintf (ferr, "Trying file %s, t = %f\n",name,time_restart);
        if( !restore(file = name)){
            fprintf(stdout,"error on restore :(\n");
            exit(1);
        }
        else{
            run();
        }
    }
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

//event testThin(t += t_out; t <= T_END){
//    foreach(){
//        phii[] = 2*f[] - 1;
//        sign[] = 7;
//    }
//    int l_sign = sign_lev;
//    for(int ilev = depth() - 1; ilev >= l_sign; ilev--){
//        foreach_level(ilev){
//            if(is_refined(cell))
//            restriction_average(point,phii);
//        }
//    }
//    compute_signature_neigh_level(f,phii,sign,l_sign);
//    int countnow = 0;
//    //foreach_level(l_sign){
//    //    //set values for each level 
//    //    int val = sign[];
//    //    if(val == -1 || val == 2){
//    //        sign[] = -1;
//    //    }
//    //    else if (val == 0){
//    //        sign[] = 0;
//    //    }
//    //    else{
//    //        sign[] = 1;
//    //        //here we add to our count
//    //        countnow++;
//    //    }
//    //}
//    if(countnow > 0){
//        fprintf(stdout,"THIN DETECTION\n");
//        hasSkeleton = true;
//    }
//    for(int ilev = l_sign; ilev < depth(); ilev++){
//        foreach_level(ilev){
//            sign.prolongation = phii.prolongation = refine_injection;
//            if(is_refined(cell)){
//                sign.prolongation(point,sign);
//                phii.prolongation(point,phii);
//            }
//        }
//    }
//}
event runSkeleton(t += t_out; t<= T_END){
    fprintf(stdout,"Skeleton at %f\n",t);
    
    //output mpi stats
    foreach(){
        int rank,size;
        rank = pid();
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        fprintf(stdout,"(%d,%d) [%f,%f] = %f\n",rank,size,x,y,f[]);
    }








    //setup
    //char sname[80];
    //sprintf (sname, "smoothskeleton-%5.3f.dat", t);
    //
    //double alpha = 25 * PI / 180;/////INPUT ANGLE TO SKELETON 
    //double **skeleton;
    //struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    //int snr;int snd;
    //double calc_time = 0;
    //
    //char redname[80];
    //sprintf(redname, "reducedskeleton-%5.3f.dat", t);
    //
    //clock_t begin = clock();
    ////calc minimumdistance
    //double mindis = 10.;
    //foreach(serial){
    //    if(Delta < mindis){
    //        mindis = Delta;
    //    }
    //}
    //fprintf(stdout,"thresh:%f\n",mindis); 
    //
    //
    ////run smooth skeleton 
    //double **sinterface = output_points_2smooth(sP,&snr,&snd,t);
    //skeleton = skeletize(sinterface,&snr,&snd,sname,&mindis,false,alpha);
    //
    //clock_t end = clock();
    //calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
    //fprintf(stdout,"time took for smooth skeleton: %f\n",calc_time);
    //
    ////next we thin it out on alpha & mindis
    //double thindis = 3*mindis;
    //thinSkeleton(&skeleton,&snd,&snr,&alpha,&thindis);

    //calc_time = 0;
    //begin = clock();
    //
    ////setup spline
    //int skelen = 3;//n of splines
    //double del = L/pow(max_level,2);
    //double minbranchlength = 0.01;
    //int mxpt = 100; 
    //skeleReduce(skeleton,del,&minbranchlength,&snr,&snd,&mxpt,t,skelen);
    //
    //end = clock();
    //calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
    //fprintf(stdout,"time took for reduced skeleton: %f\n",calc_time);
    ////cleanup data
    //for(int i = 0; i < snr;i++){
    //    free(skeleton[i]);
    //}
    //free(skeleton);
    //skeleton = NULL;

    //fflush(ferr);
    //fprintf(stdout,"\n");
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
    
    squares(color = "level",min = 1, max = max_level);
    box();
    cells();
    draw_vof("f");
    sprintf(name, "Img/levels-s-%d-%5.3f.png",(int)max_level,t);
    save(file = name);
    
    clear();
    
    //squares(color = "sign",map = black_body,linear = false,max = 2, min = -1);
    //box();
    ////cells(lc = {0,255,255});
    //draw_vof("f");
    //sprintf(name, "Img/sign-s-%d-%5.3f.png",(int)max_level,t);
    //save(file = name);
    //
    //clear();
    
    fprintf(stdout,"done\n");
}

event snapshot (t += t_out; t<=T_END) {
    //Here we calculate time taken
    char name[80];
#if 1
    sprintf (name, "dumps/snapshot-s-%d-%5.3f.gfs",(int)max_level, t);
    output_gfs (file = name, t=t, list = {f,u,p});
#endif

    sprintf (name, "dumps/dump-s-%d-%5.3f",(int)max_level, t);
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

/**
Adapt mesh based on the volume fraction. 
*/

event adapt (i=10;i++) {
    adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, minlevel = 9, maxlevel = max_level);
    //adapt_wavelet2((scalar *){f,u.x,u.y},(double []){femax,uemax,uemax},(int []){max_level, max_level-2, max_level-2},3);
}
