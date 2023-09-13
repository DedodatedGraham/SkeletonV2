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
#include "sandbox/lchirco/signature.h"
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
        time_restart = atof (argv[1]);
        time_restart = time_restart/100;
    }
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
//if want no drops
//event droplets(t += t_out){
//    remove_droplets(f,3,1e-3);
//}
event initSkele(t += t_out; t<=T_END){
    //Here we want to initalize a skeleton in all 'thin' regions of the simulation 
    //in this test we will use the manifold death signature algorithm to detect where to place skeletons
    //first compute phii 
    foreach(){
        phii[] = 2*f[] - 1;
        s[] = 7;
    }
    //int computelevel = find_moments_level(f,,);
    int computelevel = 8;
    for (int ilev = depth() - 1; ilev >= computelevel; ilev--){  
        foreach_level(ilev){
            if(is_refined(cell)){
                restriction_average(point, phii);
            }
        }
    }
    compute_signature_neigh_level(f,phii,s,computelevel);
    for (int ilev = computelevel; ilev < depth(); ilev++){
        foreach_level(ilev){
            s.prolongation = phii.prolongation = refine_injection;
            if(is_refined(cell)){
                s.prolongation (point, s);
                phii.prolongation (point, phii);
            }
        }
    }
    //now s should contain our signatures, which we will parse and find regions we wish to skeletize :)    
    
}

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
    double cpu_time_used = ((double) (timernow - timerlast)) / (CLOCKS_PER_SEC * omp_get_num_threads());
    timertotal += cpu_time_used;
    fprintf(stdout,"Total Time used: %fm: %fs\n",floor(timertotal / 60.),timertotal - 60 * (floor(timertotal / 60.)));
    fprintf(stdout,"ran @ t=%f took (%fm : %fs)\n\n",t,floor(cpu_time_used / 60.),cpu_time_used - 60 * (floor(cpu_time_used / 60.)));
    timerlast = timernow;
}


int fov1 =22;
double tx1=-0.5;
double ty1=-0.5;
double vmax=9.129;
double pmax=0.5;
int image_width=1000;
int image_height=1000;

scalar vm[],vml[];

event images(t+=t_out){
    
    foreach(){
        vm[] = (1-f[])*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
	    vml[]= f[]*sqrt(u.x[]*u.x[]+u.y[]*u.y[]);
    }
    boundary({vm,vml});
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
    squares("f", min = 0, max = 1);
    draw_vof("f");
    sprintf(name, "vof-%5.3f.png",t);
    save(file = name);
    clear();
    
    box();
    squares(color = "level", min = 6, max = max_level);
    draw_vof("f");
    sprintf(name, "levels-%5.3f.png",t);
    save(file = name);
    clear();
    
    //create plotting scalar
    //box();
    //squares(color = "plotS",map = black_body,min = 0,max  = 2);
    //draw_vof("f");
    //sprintf(name, "thin-%5.3f.png",t);
    //save(file = name);
    //clear();
    
    box();
    squares(color = "s",map = black_body,min = -1, max = 2);
    draw_vof("f");
    sprintf(name, "splt-%5.3f.png",t);
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
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, maxlevel = max_level);
  //adapt_wavelet2((scalar *){f,u.x,u.y},(double []){femax,uemax,uemax},(int []){max_level, max_level-2, max_level-2},3);
#endif 
}

