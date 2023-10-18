/** aerobreakup of a drop 
 */
#include <time.h>
//#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
#include "two-phase.h"
//#include "two-phase-s.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "skele/skeleton.h"
#include "colormap.h"
//#include "sandbox/lchirco/signature.h"

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1

//draw
double max_level = 9;
double L = 12.;
double t_out = 0.001;       
//double T_END = 5.0;
double T_END = 0.8;
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
    timerstart = clock();
    timerlast = clock();
    printf("checking %f\n",time_restart);
    size (L);
    origin(0,0);
    init_grid (128);

    /**
    The dimensionless parameters can then be computed based on 
    rho_l, sigma, and D. */
    rho1 = 1., rho2 = rho_V; 
    mu1 = mu_L, mu2 = mu_V;
    f.sigma = 1.;
    sprintf(restartname, "dump");
 
    
    //TOLERANCE = 1.e-5; 
    //NITERMIN = 5;
    run();
}

/**
The initial drop is spherical. */
double x0 = 2.; double Rd = 0.5; 


event init (t = 0){
    char dumpname[80];
    sprintf(dumpname,"dumps/dump-%5.3f",0.650);
    if (!restore (file = dumpname)) {
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

bool hasSkeleton = false;
double tthin = 10.10;
//event testThin(t += t_out; t <= T_END){
//    if(t > tthin){        

//        foreach(){
//            phii[] = 2*f[] - 1;
//            sign[] = 7;
//        }
//        int l_sign = sign_lev;
//        for(int ilev = depth() - 1; ilev >= l_sign; ilev--){
//            foreach_level(ilev){
//                if(is_refined(cell))
//                restriction_average(point,phii);
//            }
//        }
//        compute_signature_neigh_level(f,phii,sign,l_sign);
//        int countnow;
//        foreach_level(l_sign,reduction(+:countnow)){
//            //set values for each level 
//            int val = sign[];
//            if(val == 1){
//                //here we add to our count
//                countnow++;
//            }
//        }
//        if(countnow > 0){
//            fprintf(stdout,"THIN DETECTION\n");
//            hasSkeleton = true;
//        }
//        for(int ilev = l_sign; ilev < depth(); ilev++){
//            foreach_level(ilev){
//                sign.prolongation = phii.prolongation = refine_injection;
//                if(is_refined(cell)){
//                    sign.prolongation(point,sign);
//                    phii.prolongation(point,phii);
//                }
//            }
//        }
//
//        //example build of the skelevof
//        foreach(){
//            if(countnow == 0){
//                skelevof[] = 0;
//            }
//            else{
//                //Here we find cells that will be tranisitional
//                if(sign[] == 1){
//                    active_PID[cell.pid] = 1;
//                    //thin region so we set skelevof to 2
//                    skelevof[] = 2;
//                }
//                else{
//                    skelevof[] = 0;
//                }
//            }
//        }
//
//        //finally if we are in a '0' cell and there is a 2 nearby, we will set to 1;
//        foreach(){
//            //we have the point we are looking at
//            if(skelevof[] == 0){
//                bool detect = false;
//                foreach_neighbor(1){
//                    if(skelevof[] == 2){
//                        detect = true;
//                    }
//                }
//                if(detect){
//                    skelevof[] = 1;
//                }
//            }
//        }
//    }
//}

//char plotname[80];
//
//event plots (t += t_out; t<=T_END){
//    fprintf(stdout,"plotting\n");
//    double tx = 0.;
//    double vol = 0;
//    foreach(reduction(+:vol) reduction(+:tx)){
//        double dvf =f[]*dv();
//        if(f[] > 1e-6){
//            vol+= dvf;
//            tx += dvf*x;
//        }
//    }
//    tx /= vol;
//    fprintf(stdout,"%f\n",tx);
//    view(fov=8,tx=-(tx/L),ty=-0.5,sx=1,sy=1,camera="front",width=1000,height=1000);
//    
//
//    clear();
//    
//    draw_vof("f");
//    squares(color = "level",min = 1, max = max_level);
//    box();
//    cells();
//    sprintf(plotname, "Img/levels-s-%d-%5.3f.png",(int)max_level,t);
//    save(file = plotname);
//    
//    clear();
//    
//    draw_vof("f");
//    squares(color = "skelevof",map = black_body,min = 0, max = 2);//0 for none, 1 for transition, 2 for skeleton
//    box();
//    cells();
//    sprintf(plotname, "Img/svof-s-%d-%5.3f.png",(int)max_level,t);
//    save(file = plotname);
//    
//    clear();
//    
//    //squares(color = "sign",map = black_body,linear = false,max = 2, min = -1);
//    //box();
//    ////cells(lc = {0,255,255});
//    //draw_vof("f");
//    //sprintf(plotname, "Img/sign-s-%d-%5.3f.png",(int)max_level,t);
//    //save(file = plotname);
//    //
//    //clear();
//    
//    fprintf(stdout,"done\n");
//}

/**
Adapt mesh based on the volume fraction. 
*/
event adapt (i=10;i++) {
    adapt_wavelet ({f,u.x,u.y}, (double[]){femax,femax,uemax,uemax}, minlevel = 4, maxlevel = max_level);
    //adapt_wavelet2((scalar *){f,u.x,u.y},(double []){femax,uemax,uemax},(int []){max_level, max_level-2, max_level-2},3);
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
    fprintf(stdout,"Total Time used: %fm: %fs\n",floor(timertotal / 60.),timertotal - 60 * (floor(timertotal / 60.)));
    fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
    int ncells = 0;
    foreach(reduction(+:ncells)){
        ncells++;
    }
    fprintf(stdout,"Ncells --> %d\n",ncells);
    timerlast = timernow;
}



#if _MPI
//TEST FUNCTIONS
event testMPI(t += t_out){
    int dim = 2;
    printf("start-%f @ %d\n",t,pid());
    double alpha = 25 * PI / 180;/////INPUT ANGLE TO SKELETON 
    
    //search for interface cells
    foreach(){
        if(!active_PID[cell.pid]){
            if(f[] > 1e-6 && f[] < 1.-1e-6){
                active_PID[cell.pid] = 1;
            }
        }
    }
    

    //Calc skeleton will scan skeleton vof field and determine if new skeleton needs to be placed
    double **skeleton;
    int length;
    calcSkeletonMPI(f,&alpha,&dim,max_level,L,t,&skeleton,&length,active_PID);
    for(int q = 0; q < length; q++){
        free(skeleton[q]);
    }
    free(skeleton);
    printf("end-%f @ %d\n",t,pid());
}
#endif
