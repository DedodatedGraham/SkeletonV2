#include "navier-stokes/centered.h"
#include "fractions.h"
#include "curvature.h"
#include <time.h>
#include "view.h"


int maxlevel = 10;
int initlevel = 10;
double Re = 160;
double D = 1.;
double u0 = 1.;

double t_out = 0.1;
double t_end = 0.1;
//double t_end = 10.;

face vector sf[];

u.n[left] = dirichlet ((u0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

scalar Fish[];

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    Re = atof (argv[2]);
  L0 = 6;
  origin (-L0/2, -L0/2, -L0/2);
  const face vector muc[] = {0.00078125,0.00078125,0.00078125};
  mu = muc;
  run();
}

coord oloc = {0.,0.,0.};
coord vc = {0.,0.,0.}; // the velocity of the cylinder
double rr = 1.1019*sq(.12);
double chord = 1.1;

double trailing(double ax, double ay, double az, double t,double curve){
  if(ax < 0)return 1.;
  double thickness = 5*0.12*((0.2969*(sqrt(ax)))-(0.1260*ax)-(0.3516*(sq(ax)))+(0.2843*(cube(ax)))-(0.1036*(pow(ax, 4.))));
  double top = -ay-thickness-2*curve,bot = ay-thickness+2*curve;
  return max(top,bot);
}
double trailingcap(double ax, double ay, double az,double sz, double t,double curve){
  if(ax < 0)return 1.;
  double thickness = 5*0.12*((0.2969*(sqrt(ax)))-(0.1260*ax)-(0.3516*(sq(ax)))+(0.2843*(cube(ax)))-(0.1036*(pow(ax, 4.))));
  double top = -ay-thickness-2*curve,bot = ay-thickness+2*curve;
  return(max(max(top,bot),sq(curve)+sq(ay+2*curve)+sq(az-sz)-sq(thickness)));
}
double fish(double x, double y, double z,double t){
  double tb = 0.5,bb = -0.5;//bounds
  double ax = x-oloc.x,ay = y-oloc.y, az = z-oloc.z;//location to airfoil center
  double curve = 0.05*sin ((5.15*ax) - (6.283)*t);
  //geometry
  double leadingedge = sq(ax-rr)+sq(ay+2*curve)-sq(rr);//leading edge cylinder
  double bp = z + oloc.z - tb,bm = z + oloc.z - bb;
  double b = max(bp,-bm);//create bounds
  double trailingedge = trailing(ax,ay,az,t,curve);//profile of airfoil
  leadingedge = max(b,leadingedge);
  trailingedge = max(b,trailingedge);
  leadingedge = min(leadingedge,sq(ax-rr)+sq(ay+2*curve)+sq(az-tb)-sq(rr));
  leadingedge = min(leadingedge,sq(ax-rr)+sq(ay+2*curve)+sq(az-bb)-sq(rr));
  trailingedge = min(trailingedge,trailingcap(ax,ay,az,tb,t,curve));
  trailingedge = min(trailingedge,trailingcap(ax,ay,az,bb,t,curve));
  return min(leadingedge,trailingedge);//interection of leading and trailing edges
}

event init (t = 0) {
  refine(sq(x)+sq(y)+sq(z) < sq(2) && level < initlevel);
  fraction(Fish,fish(x,y,z,t));
  foreach()
    u.x[] = u0*Fish[];
}

event moving (i++) {
  fraction(Fish,fish(x,y,z,t));
  foreach()
    foreach_dimension()
      u.x[] = Fish[]*vc.x + (1. - Fish[])*u.x[];
  boundary ((scalar *){u});
}

event dump(t += t_out){
  char dumpout[80];
  sprintf(dumpout,"dumps/dump-%1.1f",t);
  dump(file = dumpout);
}

//clock_t cmark;
event outData(t += t_out; t < t_end){
  printf("calc @ %f\n",t);
}

event tout(t+=t_out){
  char outpath[80];
  sprintf(outpath,"dat/dat-%1.1f.txt",t);
  FILE *fp = fopen(outpath,"w");
  foreach(){
    if(Fish[] > 0. && Fish[] < 1.){
      fprintf(fp,"%f %f %f\n",x,y,z);
    }
  }
  fflush(fp);
  fclose(fp);
}

char name[80];
double fov1 =5;
double tx1=0.0;
double ty1=-0.5;
double vmax=180.;
double vmaxinfc=36.;
double pmax=0.5e6;
double l2max=1e-3;
int image_width=1800;
int image_height=1200;
event plot(t += t_out){
  for(int q = 0.; q < 100; q ++){
    printf("draw @ %f\n",(double)q/10.);
    fraction(Fish,fish(x,y,z,(double)q/10.));
    view (fov = fov1, tx=tx1*0.8, ty=0.00001, quat = {0,0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("Fish", min = 0, max = vmaxinfc);
    sprintf (name, "Img/vof-vel-frontview_zoom_%03d.png",q);
    save(file=name);
    
    view (fov = fov1, tx=tx1*0.8, ty=0.00001, quat = {0,-0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("Fish", min = 0, max = vmaxinfc);
    sprintf (name, "Img/vof-vel-backview_zoom_%03d.png",q);
    save(file=name);
    printf("done @ %f\n",(double)q/10.);
  }
}

//event movefish(i++; t=0.){
//  foreach_dimension(){
//    if(dt < 1.0 && dt > 0.){
//      oloc.x += vc.x * dt;
//    }
//  }
//}


event adapt(i++){
  adapt_wavelet ((scalar*){u}, (double[]){3e-2,3e-2,3e-2}, minlevel = maxlevel-2, maxlevel = maxlevel);
}

