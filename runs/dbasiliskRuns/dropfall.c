/** aerobreakup of a drop 
 */
#include "axi.h"
#include "navier-stokes/centered.h"
#include "./two-phase_nofilter.h"
//#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
//#include "adapt2.h"
#include "../../src/skeleBasilisk.h"

#define LARGE 1e36

double max_level = 8;
double L = 8.;
double t_out = 0.01;       
//double t_end = 0.00;
double t_end = 7;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.2e-3;
//WE8.5
//double mul=9.21e-3;
//double mug=1.580e-4;
//WE35
double mul=4.54e-3;
double mug=7.81e-5;
double u0 = 0;        //free stream velocity
double h   = 0.2;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
double maxruntime = 60;
//boundary for falling drop
u.t[right] = dirichlet(0);
u.t[left] = dirichlet(0);

p[right] = dirichlet(0);
p[left]  = neumann(0);

pf[right] = dirichlet(0);
pf[left]  = neumann(0);

uf.n[bottom] = 0.;
uf.n[top] = 0.;




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
 
  run();
}

/**
The initial drop is spherical. */
double x0 = 7.; double Rd = 0.5; 

event init (t = 0){
    if (!restore (file = "dump")) {
        refine (  sq(y-x0)+sq(x) < sq(1.2*Rd) && sq(y-x0)+sq(x) > sq(0.8*Rd) && level < max_level);
        fraction (f, -sq(y-x0)-sq(x)+sq(Rd));

        /** It is important to initialize the flow field. If u.x=0, Poisson
        solver crahses. If u.x=u0*(1-f[]), the interface velocity is too large
        and will induce erroneous deformation at early time and wrong pressure 
        field. As along as we use u.x=a*u0*(1-f[]), with a<0.01, then the results 
        look normal and not sensitive to a. */
        //foreach() {
        //    u.y[]=1.e-7*u0*(1.-f[]);
        //}
        //boundary ({f,u.y});
    }
}
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] -= 0.98;
  boundary({p});
}

int slevel = 0.;
double mindis = 0.;
event skeleton(t+=t_out){
    //First find min grid distance    
    if(slevel == 0 || slevel < max_level){
        foreach(){
            if(f[] > 0. && f[] < 1.)
            if(Delta < mindis || mindis == 0.){
                mindis = Delta;
                slevel = max_level;
            }
        }
    }
    fprintf(stdout,"Skeleton at %f\n",t);
    //fprintf(stdout,"mindis= %f\n",mindis);
    //fprintf(stdout,"slevel= %d\n",slevel);
    char sname[80];
    sprintf (sname, "skeleton-%5.3f.dat", t);
    struct OutputXYNorm sP; sP.c = f; sP.level = max_level;
    int snr;int snd;
    double **sinterface = output_points_xynorm(sP,&snr,&snd);
    //fprintf(stdout,"%d %d\n",snr,snd);
    skeletize(sinterface,&snr,&snd,sname,&mindis);
    //free(sinterface);
    fflush(ferr);
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
	    double abs = sqrt(pow(n.x,2)+pow(n.y,2));
        double tx = n.x/abs;
        double ty = n.y/abs;
	    fprintf(p.fp, "%g %g %g %g\n",x+Delta*pc.x, y+Delta*pc.y, tx,ty);
	    fprintf(p.fp, "%g %g %g %g\n",-(x+Delta*pc.x), y+Delta*pc.y, -tx, ty);
	    //fprintf(p.fp, "%g %g %g %g\n",x+Delta*pc.x, -y-Delta*pc.y, n.x,-n.y);
	}
    }
    fflush(p.fp);
}


//#if !_MPI

// output snapshot of simulation 
event interface (t += t_out) {

    char name[80];
    sprintf (name, "infc-%5.3f.dat", t);
    FILE * fp1 = fopen (name, "w");
  
    output_points_norm(f,fp1,level=max_level);
  
    //output_facets (f,fp1);
    fflush(fp1);
    fclose(fp1);
}
//#endif 
 
event snapshot (t += t_out; t<=t_end ) {
    char name[80];
#if 1
    sprintf (name, "snapshot-%5.3f.gfs", t);
    output_gfs (file = name, t=t, list = {f,u,p});
#endif

    sprintf (name, "dump-%5.3f", t);
    p.nodump = false;
    dump (file = name); // so that we can restart
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
#endif 

}

