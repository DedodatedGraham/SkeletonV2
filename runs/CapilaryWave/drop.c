#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "view.h"
#include "tag.h"
#include "embed.h"

double max_level = 7;
double femax = 0.001;
double uemax = 0.001;

double L = 8 [1];
double t_out = 0.010 [0,1];       
double t_end = 0.3001 [0,1];    

double r = 0.5;
int NumDrop = 2;
/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.297E-03[*];
double mul=5.275E-03;
double mug=9.077E-05;
double u0 =1.964E+01;        //free stream velocity
double Re;
double hx = -2.0; 
double hy = -2.0; 
double maxruntime = 14420;

u.n[left] = dirichlet(u0);
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);

int main(int argc, char * argv[]){
  size(L);
  origin(hx,hy);
  init_grid(128);
  periodic(right);
  rho1 = 1. [*], rho2 = rhog ;
  mu1 = mul, mu2 = mug;
  f.sigma = 1.;
  TOLERANCE = 1.e-5 [*];
  run();
}

event init(t = 0.){
    refine ( level < max_level );
    for(int i = 0; i < NumDrop; i++) {
        scalar tempf[] = f[];
        fraction(tempf, -sq(x)-sq(y)+sq(r));
        foreach(){
            f[] = max(tempf[],f[]);
        }
    }
    //fraction(f,-sq(x)-sq(y)+sq(r));


    double d = HUGE; 
    foreach(){
        u.x[]=1.e-7*u0*(1.-f[]);
	    if (Delta < d) {
	    	d = Delta;
	    }
    }
    boundary({f,u.x});
}

char outpath[80];
event plot(t += t_out; t <= t_end){
     view(tx=-0.25,ty=-0.25,sx=1,sy=1,camera="front",width=1000,height=1000);
     clear();
     sprintf(outpath,"Img2/ux-%5.3f.png",t);
     draw_vof("f");
     squares(color = "u.x");
     box();
     save(file = outpath);
     clear();

     sprintf(outpath,"Img2/vof-%5.3f.png",t);
     draw_vof("f");
     squares(color = "f");
     box();
     save(file = outpath);
     clear();
}

event outputFlowfield(t += t_out; t <= t_end){
    sprintf(outpath,"datregion/dropinformation-%5.3f",t);
    FILE *fp = fopen(outpath,"w");
//    foreach(){
//        Re=(u.x[]*2*r)/mul;
//        fprintf(fp,"%f %f %f %f %f\n", f[], u.x[],u.y[],p[],Re);
//    }
    double xc = 0.;
    double yc = 0.;
    int centcount = 0;
    foreach(){
        if(f[]){
            xc += x;
            yc += y;
            centcount++;
        }
    }
    double halflength = 1.5;
    xc = xc / centcount;
    yc = yc / centcount;
    coord pc;
    coord bounds[2] = {{ xc - halflength, yc - halflength }, {xc + halflength, yc + halflength}};
    coord nn = {128, 128};
    foreach_region(pc, bounds, nn) {
       Re=(u0*2*r)/mul;
       fprintf(fp,"%f %f %f %f %f\n", f[], u.x[],u.y[],p[],Re);
    }
    fflush(fp);
    fclose(fp);
}

//event adapt (i=10;i++) {
//  adapt_wavelet ({f,u}, (double[]){femax,uemax,uemax}, minlevel = 7, maxlevel = max_level);
//}

