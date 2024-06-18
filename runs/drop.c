/** aerobreakup of a drop 
 */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"
#include "view.h"
#include "axi.h"

double max_level = 10;
double L = 5 [1];
double t_out = 0.010 [0,1];       
double t_end = 0.3001 [0,1];    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.081e-3 [*];
double mul=4.232e-2 [*];
double mug=4.731e-5 [*];
double u0 =118.16;        //free stream velocity
double hx = -1.5 ; 
double hy = -2.5 ; 
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 14420;
double t_amr = 0.0; 

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
    hx        = atof (argv[6]);
  if (argc > 7)
    hy        = atof (argv[7]);
  if (argc > 8)
    rhog      = atof (argv[8]);
  if (argc > 9)
    mul       = atof (argv[9]);
  if (argc > 10)
    mug       = atof (argv[10]);
  if (argc > 11)
    femax     = atof (argv[11]);
  if (argc > 12)
    uemax     = atof (argv[12]);
  if (argc > 13)
    maxruntime= atof (argv[13]);
  if (argc > 14)
    t_amr     = atof (argv[14]);
  #if _OscilDrop
  if (argc > 13)
    A2        = atof (argv[13]);
  if (argc > 14)
    A3        = atof (argv[14]);
  #endif

  size (L);
  init_grid (128);

  /**
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */

  rho1 = 1. [*], rho2 = rhog ; 
  mu1 = mul, mu2 = mug;
  f.sigma = 1.;
  TOLERANCE = 1.e-5 [*]; 
  //VISCOUS_TOL = 4.e-3; 
  NITERMIN = 2 [*];
  NITERMAX = 50 [*];
  run();
}

/**
The initial drop is spherical. */
event init (t = 0)
{
  if (!restore (file = "dumpNONE")) {
    double Rd = 0.5;
    int init_max_level = max_level; 
    refine (  sq(y)+sq(x - 4*Rd)<sq(0.6) && sq(y)+sq(x - 4*Rd)>sq(0.4) && level < init_max_level );
    fraction (f, -sq(y)-sq(x - 4*Rd)+sq(Rd));
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

event removespurious(i++){
    remove_droplets(f,3,0);
    remove_droplets(f,3,0,true);
}

#include "skele/skeleInclude.h"
event makeSkelePost(t += t_out){
  double alpha=0.7;
  double **skeleton = NULL;//skeleton container
  int skelelength = 0;
#if _MPI
  calcSkeletonMPI(f,&alpha,max_level,L0,t,&skeleton,&skelelength,0);
#else
  //make interface compare 
  //non mpi only
  calcSkeleton(f,&alpha,max_level,L0,t,&skeleton,&skelelength,0);
  scalar skelef[] = f[];
  foreach(){
    skelef[] = 0.;
  }
  for(int q = 0; q < skelelength; q++){
    //load skelef with a max fraction
    double *spt = skeleton[q];
    scalar news[];
    fraction(news,-sq(x-spt[0])-sq(y-spt[1])+sq(spt[2]));
    foreach(){
      skelef[] = max(skelef[],news[]);
    }
  } 
  char outsurface[80];
  sprintf(outsurface,"dat/outsurface.txt");
  FILE *fouts = fopen(outsurface,"a");
  double farea = interface_area(f);
  double sarea = interface_area(skelef);
  double fvol = 0.;
  double svol = 0.;
  foreach(){
      fvol = fvol + f[];
      svol = svol + skelef[];
  }
  fprintf(fouts,"%5.3f %f %f %f %f\n",t,farea,sarea,fvol,svol);
  fflush(fouts);
  fclose(fouts);
#endif
  for(int q = 0; q < skelelength; q++)free(skeleton[q]);
  free(skeleton);
}

event snapshot (t += t_out; t<=t_end ) {
  printf("calc @ t=%5.3f, n = %ld\n",t,grid->tn);
  char name[80];
  sprintf (name, "dumps/dump-%06.4f", t);
  p.nodump = false;
  dump (file = name); // so that we can restart
}

/**
Adapt mesh based on the volume fraction. 
*/
event adapt (t=t_amr;i++) {
  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, minlevel = 7, maxlevel = max_level);
}

#if TRACE > 1
event profiling (i += 20) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif

event logoutput(t += t_out){
    double tx=0.,ty=0.;
    int n = 0;
    foreach(reduction(+:tx) reduction(+:ty) reduction(+:n)){
        if(f[] > 1e-6 && f[] < 1-1e-6){
            tx = tx + x;
            ty = ty + y;
            n++;
        }
    }
    tx = tx / n;
    ty = ty / n;
#if _MPI
    if(pid() == 0)
#endif
    printf("drop position [%f,%f] from %d cells\n",tx,ty,n);
}
