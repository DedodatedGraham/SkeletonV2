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

//#define T_ADAPT_OFF 1
#define LARGE 1e36
//#define EVAP_OFF 1
double max_level = 10;
double L = 12.;
double t_out = 0.01;       
double T_END = 0.45;    
//double T_END = 0.1;
//double T_END = 0.01;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rho_L=1;
double rho_V=1.297e-3;
double mu_L=5.275e-3;
double mu_V=9.077e-5;
double u0 = 87.8;        //free stream velocity
//double u0 = 107.5;
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 60;

double time_restart = 0.;
//timer
clock_t timerstart;
clock_t timernow;
clock_t timerlast;

//#include "01_vaporization/evap_include.h"

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);


int main(int argc, char * argv[])
{

  if (argc > 1)
    time_restart = atof (argv[1]);
  timerstart = clock();
  timerlast = clock();
  char name[80];
  if(time_restart == 0.0){
    size (L);
    origin(0,0);
    init_grid (64);
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
void output_skeleinterface(char name[80],double **list,int length){
    FILE *fp = fopen(name,"w");
    for(int i = 0; i < length; i++){
        //fprintf(stdout,"\n(%d / %d)\n",i,length);
        //fprintf(stdout,"Outputting [%f,%f]\n",list[i][0],list[i][1]);
        //fprintf(stdout,"[%f,%f]\n",list[i][2],list[i][3]);
        fprintf(fp,"%f %f %f %f\n",list[i][0],list[i][1],list[i][2],list[i][3]);
    }
    fflush(fp);
    fclose(fp);
}

//#include "01_vaporization/adapt_evap.h"
//Function For Obtaining Skeleton
//double mindis = 0.0;
//int slevel = 0.;
//event skeleton(t+=t_out){
//    #pragma omp master
//    {
//        fprintf(stdout,"Talking from %d\n",pid());
//        if(slevel == 0 || slevel < max_level){
//            foreach(serial){
//                if(f[] > 1e-6 && f[] < 1-1e-6)
//                if(Delta < mindis || mindis == 0.){
//                    mindis = Delta;
//                    slevel = max_level;
//                }
//            }
//        }
//        fprintf(stdout,"thresh:%f\n",mindis); 
//        fprintf(stdout,"Skeleton at %f\n",t);
//        
//        //setup smooth
//        struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
//        int snr;int snd;
//        
//        //run smooth
//        char sname[80];
//        sprintf (sname, "smoothskeleton-%5.3f.dat", t);
//        double calc_time = 0;
//        clock_t begin = clock();
//         
//        double **sinterface = output_points_2smooth(sP,&snr,&snd,t);
//        skeletize(sinterface,&snr,&snd,sname,&mindis,false);
//        
//        clock_t end = clock();
//        calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
//        fprintf(stdout,"time took for smooth skeleton: %f\n",calc_time);
//        
//        //run VOF
//        char vsname[80];
//        sprintf (vsname, "vofskeleton-%5.3f.dat", t);
//        calc_time = 0;
//        begin = clock();
//        
//        double **vsinterface = output_points_xynorm(sP,&snr,&snd);
//        skeletize(vsinterface,&snr,&snd,vsname,&mindis,true);
//        
//        end = clock();
//        calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
//        fprintf(stdout,"time took for vof skeleton: %f\n",calc_time);
//        
//        //Finally run Level Set Method
//        //char lsname[80];
//        //sprintf (lsname, "LSskeleton-%5.3f.dat", t);
//        //calc_time = 0;
//        //begin = clock();
//        //
//        ////reinit_LS(f);
//        //end = clock();
//        //calc_time = calc_time + (double)(end-begin)/CLOCKS_PER_SEC;// this is the time required for skeleton  
//        //fprintf(stdout,"time took for ls skeleton: %f\n",calc_time);
//        
//        fflush(ferr);
//    
//    
//    }
//}
//
////Function for extracting the interfacial points: cut-surface center
//struct OutputPoints{
//    scalar c;
//    FILE *fp;
//    face vector s;
//    int level;
//};
//
//
//void output_points_norm(struct OutputPoints p){
//    scalar c = p.c;
//    restriction({c});
//    face vector s = p.s;
//    if(!p.fp) p.fp = stdout;
//    if(!s.x.i) s.x.i = -1;
//
//    //fprintf(p.fp, "#1:x y\n");
//    int j = 0;
//    foreach(serial){
//        if(c[] > 1e-7 && c[] < 1.-1e-7){
//	    coord n = facet_normal(point, c, s);
//	    double alpha = plane_alpha(c[],n);
//	    coord pc;
//	    double area = plane_area_center(n, alpha, &pc);
//	    if(area==0){
//	        fprintf(stdout,"Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
//	    }
//	    double abs = sqrt(pow(n.x,2)+pow(n.y,2));
//        double tx = n.x/abs;
//        double ty = n.y/abs;
//	    fprintf(p.fp, "%g %g %g %g\n",x+Delta*pc.x, y+Delta*pc.y, tx, ty);
//	    //fprintf(p.fp, "%g %g %g %g\n",x+Delta*pc.x, -(y+Delta*pc.y), tx, -ty);
//	    //fprintf(p.fp, "%g %g %g %g\n",x+Delta*pc.x, -y-Delta*pc.y, n.x,-n.y);
//        j++;
//	    }
//    }
//    fflush(p.fp);
//}


//#if !_MPI
// output snapshot of simulation 
//event interface (t += t_out) {
//
//    char name[80];
//    sprintf (name, "infc-%5.3f.dat", t);
//    FILE * fp1 = fopen (name, "w");
//  
//    output_points_norm(f,fp1,level=max_level);
//  
//    //output_facets (f,fp1);
//    fflush(fp1);
//    fclose(fp1);
//}

//#endif 
event snapshot (t += t_out; t<=T_END) {
    //Here we calculate time taken
    char name[80];
#if 1
    sprintf (name, "snapshot-%5.3f.gfs", t);
    output_gfs (file = name, t=t, list = {f,u,p});
#endif

    sprintf (name, "dump-%5.3f", t);
    p.nodump = false;
    dump (file = name); // so that we can restart

    //finally output timer and results as this takes bulk of process
    timernow = clock();
    double cpu_time_used = ((double) (timernow - timerlast)) / (CLOCKS_PER_SEC * 14);
    fprintf(stdout,"%f\n",cpu_time_used);
    fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
    timerlast = timernow;
}


int fov1 =22;
double tx1=-0.5;
double ty1=-0.5;
double vmax=9.129;
double pmax=0.5;
int image_width=600;
int image_height=600;

scalar vm[],vml[];

event images(t+=t_out){
    
    foreach(){
        vm[] = (1-f[])*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
	vml[]= f[]*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
    }boundary({vm,vml});
    scalar Omega[];
    vorticity(u,Omega);
    clear();
    char name[80];
    view(fov = fov1, tx = tx1, ty = ty1, width = image_width, height =image_height);
    box();
    cells();
    draw_vof("f");
    squares("vm",min = 0, max = vmax);
    sprintf(name,"vof_velmag-%5.3f.png",t);
    save(file = name);
    clear();
	
    view(fov = fov1, tx = tx1, ty = ty1, width = image_width, height =image_height);
    //cells();
    box();
    draw_vof("f");
    squares("f", min = 0, max = 1);
    sprintf(name, "vof-%5.3f.png",t);
    save(file = name);
    clear();
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
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, minlevel = 3, maxlevel = max_level);
  //adapt_wavelet2((scalar *){f,u.x,u.y},(double []){femax,uemax,uemax},(int []){max_level, max_level-2, max_level-2},3);
#endif 
}

