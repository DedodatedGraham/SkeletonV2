#include "navier-stokes/centered.h"
#include "fractions.h"
#include "curvature.h"
#include <time.h>
#include "view.h"

int maxlevel = 10;
double Re = 160;
double D = 1.;
double u0 = 1.;

double t_out = 0.01;
double t_end = 19.9;

face vector sf[];


u.n[left] = dirichlet ((u0));
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

int main (int argc, char * argv[])
{
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    Re = atof (argv[2]);
  L0 = 10;
  origin (-L0/2, -L0/2);
  const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  run();
}

scalar Fish[];
//represents a moving fish
coord oloc = {0.,0.};
double rr = 1.1019*sq(.12);

double trailing(double ax, double ay, double t,double curve){
  if(ax < 0)return 1.;
  double thickness = 5*0.12*((0.2969*(sqrt(ax)))-(0.1260*ax)-(0.3516*(sq(ax)))+(0.2843*(cube(ax)))-(0.1036*(pow(ax, 4.))));
  double top = -ay-thickness-2*curve,bot = ay-thickness+2*curve;
  return max(top,bot);
}

double fish(double x, double y,double t){
  double ax = x-oloc.x,ay = y-oloc.y;//location to airfoil center
  double curve = 0.05*sin ((5.15*ax) - (6.283)*t);
  //geometry
  double leadingedge = sq(ax-rr)+sq(ay+2*curve)-sq(rr);//leading edge
  double trailingedge = trailing(ax,ay,t,curve);//profile of airfoil
  return min(leadingedge,trailingedge);//interection of leading and trailing edges
}

event init (t = 0.) {
  //mask(y > 0.5 ? top: y < -0.5 ? bottom : none);
  refine(level < maxlevel);
  fraction(Fish,fish(x,y,t));
  foreach()
    u.x[] = u0*(Fish[]);
}

event moving (i++) {
  fraction(Fish,-fish(x,y,t));
  //coord vc = {u0,0.,0.};
  foreach()
    foreach_dimension()
      u.x[] = (1. - Fish[])*u.x[];
  boundary ((scalar *){u});
}

char outpath[80];
event plot (t += t_out; t < t_end){
  
  scalar omega[];
  vorticity (u, omega);
  


  view(tx=0,ty=0,sx=1,sy=1,camera="front",width=1000,height=1000);
  clear();
  sprintf(outpath,"Img/omega-%5.3f.png",t);
  draw_vof("Fish");
  squares(color = "omega",min = -12, max = 12);
  box();
  save(file = outpath);
  
  clear();
  sprintf(outpath,"Img/ux-%5.3f.png",t);
  draw_vof("Fish");
  squares(color = "u.x");
  box();
  save(file = outpath);
  
  clear();
  sprintf(outpath,"Img/levelint-%5.3f.png",t);
  draw_vof("Fish");
  squares(color = "level",min = 1, max = maxlevel);
  box();
  save(file = outpath);
  
  clear();
}

event outData(i++){
  printf("@ %f / %f \n",t,t_end);
  //Output important info for plotting
  scalar omega[];
  vorticity (u, omega);
  sprintf(outpath,"dat/scalarplot-%5.3f.txt",(double)i/1000.);
  FILE *fp = fopen(outpath,"w");
  foreach(){
    if(Fish[] < 1.){
      fprintf(fp,"%f %f %f %f %f %f %f\n",x,y,Delta,u.x[],u.y[],p[],omega[] > 12. ? 12. : omega[] < -12. ? -12. : omega[]);
    }
  }
  fflush(fp);
  fclose(fp);
  //output geometrical features
  sprintf(outpath,"dat/geometry-%5.3f.txt",(double)i/1000.);
  FILE *fp1 = fopen(outpath,"w");
  foreach(serial){
    if(Fish[] > 0. && Fish[] < 1.){
      //compute local drag force acting
      //first dynamic pressure
      fprintf(fp1,"%f %f\n",x,y);
    }
  }
  fflush(fp1);
  fclose(fp1);
}

#include "skele/skeleInclude.h"
event makeSkeletize(i++){
  double alpha = 20 * PI / 180;/////INPUT ANGLE TO SKELETON 
  printf("%f\n",alpha);
  double **skeleton = NULL;//skeleton container
  int skelelength = 0;
  calcSkeleton(Fish,&alpha,maxlevel,L0,(double)i/1000.,&skeleton,&skelelength,0);
  for(int q = 0; q < skelelength; q++)free(skeleton[q]);
  free(skeleton);
}

//event adapt(i++){
//  adapt_wavelet ((scalar*){u}, (double[]){3e-2,3e-2}, minlevel = maxlevel-5, maxlevel = maxlevel);
//}
