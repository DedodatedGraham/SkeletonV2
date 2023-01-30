/** aerobreakup of a drop 
 */
#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
//#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
//#include "adapt2.h"
//#include "legendre.h"
#include "../src/skeletize.h"

#define LARGE 1e36

double max_level = 8;
double L = 4.;
double t_out = 0.01;       
double t_end = 4;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.2e-3;
double mul=2.995e-3;
double mug=5.390e-5;
double u0 = 16.71;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 60;

u.n[left] = dirichlet(u0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);


int main(int argc, char * argv[])
{

  if (argc > 1)
    max_level = atoi (argv[1]);
  if (argc > 2)
    L         = atof (argv[2]);
  if (argc > 3)
    u0        = atof (argv[3]);
  if (argc > 4)
    t_out     = atof (argv[4]);
  if (argc > 5)
    t_end     = atof (argv[5]);
  if (argc > 6)
    rhog      = atof (argv[6]);
  if (argc > 7)
    mul       = atof (argv[7]);
  if (argc > 8)
    mug       = atof (argv[8]);
  if (argc > 9)
    femax     = atof (argv[9]);
  if (argc > 10)
    uemax     = atof (argv[10]);
  if (argc > 11)
    maxruntime = atof (argv[11]);

  size (L);
  origin(0,0);
  init_grid (64);

  /**
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */

  rho1 = 1., rho2 = rhog; 
  mu1 = mul, mu2 = mug;
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
        refine (  sq(y)+sq(x-x0) < sq(1.2*Rd) && sq(y)+sq(x-x0) > sq(0.8*Rd) && level < max_level);
        fraction (f, -sq(y)-sq(x-x0)+sq(Rd));

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

//Function for extracting the interfacial points: cut-surface center
struct OutputPoints{
    scalar c;
    FILE *fp;
    face vector s;
    int level;
};

void output_points_norm(struct OutputPoints p){
    scalar c = p.c;
    restriction({c});
    face vector s = p.s;
    if(!p.fp) p.fp = stdout;
    if(!s.x.i) s.x.i = -1;

    //fprintf(p.fp, "#1:x y\n");
    foreach_level_or_leaf(p.level){
        if(c[] > 1e-7 && c[] < 1.-1e-7){
	    coord n = facet_normal(point, c, s);
	    double alpha = plane_alpha(c[],n);
	    coord pc;
	    double area = plane_area_center(n, alpha, &pc);
	    if(area==0){
	        fprintf(stdout,"Area=Null\n");// This statement is just to make some use of the area info. otherwise compiler throws warning!!
	    }
	    fprintf(p.fp, "%g %g\n",x+Delta*pc.x, y+Delta*pc.y);
	}
    }
    fflush(p.fp);
}


 

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

//#if _UNREFINE
//  unrefine ( x>0.5*L && level >= max_level-2 );
//#endif

}

