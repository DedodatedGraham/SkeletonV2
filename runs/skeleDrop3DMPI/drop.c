/** aerobreakup of a drop 
 */

#include <time.h>
#include "navier-stokes/centered.h"
#include "two-phase-s.h"

#include "navier-stokes/conserving.h"
#include "tension.h"

#include "curvature.h"
#include "view.h"
#include "tag.h"
#include "skele/skeleton.h"
#include "colormap.h"

#define LARGE 1e36

double max_level = 9;
double L = 5;
double t_out = 0.010;       
double t_end = 0.5;
//double t_end = 0.320;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.081e-3;
double mul=4.232e-2;
double mug=4.731e-5;
double u0 =118.16;        //free stream velocity
double h = 2.5; 
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 14420;
double t_amr = 0.0; 

#if 0 
double gravity = 9.8;     //gravity acceleration (m/s2)
double gx;  //dimensionless gravity in x direction  
#endif


#if 0
u.n[left]  = dirichlet(u0);
u.t[left]  = dirichlet(0);
p[left]    = neumann(0);

u.n[top]  = dirichlet(0);
u.t[top]  = neumann(0);
p[top]    = neumann(0);

u.n[right]  = neumann(0);
u.t[right]  = neumann(0);
p[right]    = dirichlet(0);
#else 
u.n[left] = dirichlet(u0);

u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);
#endif 

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
    h        = atof (argv[6]);
  if (argc > 7)
    rhog      = atof (argv[7]);
  if (argc > 8)
    mul       = atof (argv[8]);
  if (argc > 9)
    mug       = atof (argv[9]);
  if (argc > 10)
    femax     = atof (argv[10]);
  if (argc > 11)
    uemax     = atof (argv[11]);
  if (argc > 12)
    maxruntime= atof (argv[12]);
  if (argc > 13)
    t_amr     = atof (argv[13]);
  #if _OscilDrop
  if (argc > 14)
    A2        = atof (argv[14]);
  if (argc > 15)
    A3        = atof (argv[15]);
  #endif

  size (L);
  origin (0, -L/2, -L/2.);
  init_grid (64);

  /**
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */

  rho1 = 1., rho2 = rhog; 
  mu1 = mul, mu2 = mug;
  f.sigma = 1.;
 
#if 0 
  gx  = gravity*1100.*sq(2.7e-3)/0.0483;
#endif 

  TOLERANCE = 1.e-5; 
  NITERMIN = 2;
  NITERMAX = 30;

  run();
}

/**
The initial drop is spherical. */
double voidintsphere(double x,double y,double z,double r,double x0, double y0, double z0){
    double outer = -sq(x-x0)-sq(y-y0)-sq(z-z0)+sq(r);
    double inner = -sq(x-x0)-sq(y-y0)-sq(z-z0)+sq((r-0.1)+ 0.07*sin(2.*atan2(sqrt(sq(y-y0)+sq(z-z0)),x-x0)));
    if(inner <= 0.){
        //not inside inner
        if(x-x0 > 0. && outer > 0.){
            //if in right half
            return outer;
        }
    }
    if(x-x0 >= -0.2 && x-x0 <= 0.2){
        //place our torus
        double tor = -sq((r-0.05) - sqrt(sq(y-y0) + sq(z-z0)))-sq(x-x0)+sq(0.05); 
        if(tor > 0.)return tor;
    }
    //if nothing else we return our inside val
    if(inner < 0.)return inner;
    return -inner;
}
event init (t = 0)
{
  if (!restore (file = "importDump/dump-0.4280")) {

    double Rd = 0.55;
    double x0 = 0.55+h;
    int init_max_level = min(12,max_level); 

    refine (  sq(y)+sq(x-x0)+sq(z)<sq(0.6) && sq(y)+sq(x-x0)+sq(z)>sq(0.4) && level < init_max_level );
    fraction (f,voidintsphere(x,y,z,Rd,x0,0.,0.));
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
  else{
      printf("loaded!\n");
  }
}

event snapshot (t += t_out; t<=t_end ) {
  char name[80];
  printf("t=%05.3f\n",t);
  sprintf (name, "dumps/snapshot-%05.3f.gfs", t);
  output_gfs (file = name, t=t, list = {f,u,p});

  sprintf (name, "dumps/dump-%05.3f", t);
  p.nodump = false;
  dump (file = name); // so that we can restart
}

/**
Adapt mesh based on the volume fraction. 
*/
event adapt (t=t_amr;i++) {
#if 0
  vector q[];
  foreach()
    foreach_dimension()
      q.x[] = u.x[]*(rho1*f[]+rho2*(1.-f[]));
  boundary ((scalar *){q});
  adapt_wavelet ({f,q}, (double[]){femax,uemax,uemax}, minlevel = 3, maxlevel = max_level);
#else
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax,uemax}, minlevel = max_level-7, maxlevel = max_level);
#endif 

#if _UNREFINE
  /** the domain is -3<x<13 if hx=-3*/
  unrefine ( (x>10. && level >= max_level-3)||(x>9. && level >= max_level-2) );
#endif
}

event runtime (i += 10) {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t/60 >= maxruntime) {
    dump (file = "dump"); // so that we can restart
    return 1; // exit
  }
}
