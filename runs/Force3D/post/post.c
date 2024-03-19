/** aerobreakup of a drop 
 */
#include <string.h>

#include "navier-stokes/centered.h"
#include "lambda2.h"
#include "view.h"
#include "iso3D.h"

#define LARGE 1e36

int i_start = 0;
int i_end   = 200;
int i_gap   = 10;
double fov1 =25;
double tx1=-0.5;
double ty1=-0.5;
double vmax=180.;
double vmaxinfc=36.;
double pmax=0.5e6;
double l2max=1e-3;
int image_width=1800;
int image_height=1200;

coord oloc = {0.,0.,0.};
coord vc = {0.,0.,0.}; // the velocity of the cylinder
double rr = 1.1019*sq(.12);
double chord = 1.1;

double trailing(double ax, double ay, double az, double t){
  double curve = 0.05*sq(ay)*sin ((5.15*ax) - (6.283)*t);
  double thickness = 5*0.12*((0.2969*(sqrt(ax)))-(0.1260*ax)-(0.3516*(sq(ax)))+(0.2843*(cube(ax)))-(0.1036*(pow(ax, 4.))));
  return -1. * (ay > curve? - ay + (curve + thickness):ay - (curve - thickness));
}
double trailingcap(double ax, double ay, double az,double sz, double t){
  double curve = 0.05*sq(ay)*sin ((5.15*ax) - (6.283)*t);
  double thickness = 5*0.12*((0.2969*(sqrt(ax)))-(0.1260*ax)-(0.3516*(sq(ax)))+(0.2843*(cube(ax)))-(0.1036*(pow(ax, 4.))));
  return sq(curve)+sq(ay)+sq(az-sz)-sq(thickness);
}
double fish(double x, double y, double z,double t){
  double tb = 0.5,bb = -0.5;//bounds
  double ax = x-oloc.x,ay = y-oloc.y, az = z-oloc.z;//location to airfoil center
  //geometry
  double leadingedge = sq(ax-rr)+sq(ay)-sq(rr);//leading edge cylinder
  double bp = z + oloc.z - tb,bm = z + oloc.z - bb;
  double b = max(bp,-bm);//create bounds
  double b1p = x + oloc.x - chord,b1m = x + oloc.x;
  double b1 = max(b1p,-b1m);//create bounds
  double trailingedge = trailing(ax,ay,az,t);//profile of airfoil
  leadingedge = max(b,leadingedge),trailingedge = max(b,trailingedge);
  leadingedge = min(leadingedge,sq(ax-rr)+sq(ay)+sq(az-tb)-sq(rr));
  leadingedge = min(leadingedge,sq(ax-rr)+sq(ay)+sq(az-bb)-sq(rr));
  trailingedge = min(trailingedge,trailingcap(ax,ay,az,tb,t));
  trailingedge = min(trailingedge,trailingcap(ax,ay,az,bb,t));
  return max(b1,min(leadingedge,trailingedge));//interection of leading and trailing edges
}

int main(int argc, char * argv[])
{
  if (argc > 1)
    i_start = atoi (argv[1]);
  if (argc > 2)
    i_end   = atoi (argv[2]);
  if (argc > 3)
    i_gap   = atoi (argv[3]);
  if (argc > 4)
    fov1    = atof (argv[4]);
  if (argc > 5)
    tx1     = atof (argv[5]);
  if (argc > 6)
    ty1     = atof (argv[6]);
  if (argc > 7)
    vmax    = atof (argv[7]);
  if (argc > 8)
    pmax    = atof (argv[8]);
  if (argc > 9)
    image_width  = atoi (argv[9]);
  if (argc > 10)
    image_height = atoi (argv[10]);
  if (argc > 11)
    l2max   = atof (argv[11]);
  if (argc > 12)
    vmaxinfc = atof (argv[12]);

  char name[80];

  int i=0; 
  for ( i = i_start; i <= i_end; i+=i_gap ) {
    //scalar l2[];
    //lambda2 (u, l2);
    
    sprintf (name, "../dumps/dump-%1.1f", (double)i/10);
    restore(file = name);
    fprintf (ferr, "Opening file %s, i=%d \n",name,i);
    scalar Fish[];
    fraction(Fish,fish(x,y,z,t));
    view (fov = fov1, tx=tx1*0.8, ty=0.00001, quat = {0,0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("Fish", min = 0, max = vmaxinfc);
    sprintf (name, "Img/vof-vel-frontview_zoom_%04d.png",i);
    save(file=name);
    
    view (fov = fov1, tx=tx1*0.8, ty=0.00001, quat = {0,-0.234,0,0.971}, width = image_width, height = image_height);
    clear();
    draw_vof ("Fish", min = 0, max = vmaxinfc);
    sprintf (name, "Img/vof-vel-backview_zoom_%04d.png",i);
    save(file=name);
    
    //view (fov = fov1*1.3, tx=tx1*2, ty=ty1, quat = {0,0.234,0,0.971}, width = image_width, height = image_height);
    //clear();
    //draw_vof ("Fish", min = 0, max = vmaxinfc);
    //isosurface ("l2", -l2max);
    //sprintf (name, "Img/vof-l2-frontview_%04d.png",i);
    //save(file=name);
  }
}

