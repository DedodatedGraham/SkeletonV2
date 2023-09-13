/** aerobreakup of a drop 
 */
#include "axi.h"
#include "navier-stokes/centered.h"
//#include "./two-phase_nofilter.h"
//#include "./two-phase_filter.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "tag.h"

#define LARGE 1e36

#if _MD
//#define SQUARES 1
#include "signature.h"
//#include "output_hdf.h"
scalar cf[], phii[], sign[], M[];
//const int max_change=1;
bool large = true;

#endif

int max_level = 10;
double L = 16.;
double t_out = 0.010;       
double t_end = 0.3001;    

/** dimensionless properties, normalized by scaling variables rhol, D, sigma
 */
double rhog=1.081e-3;
double mul=4.232e-2;
double mug=4.731e-5;
double u0 =118.16;        //free stream velocity
double h   = 2.5;          //initial gap between drop and inlet
double femax = 0.001;
double uemax = 0.001;
#if _MD
int sign_lev = 8;
#endif
int ON = 1;
#if _SlipTop
double yh = 2.5;
#endif

#if _OscilDrop     
double A2=0.0;
double A3=0.0;
double A4=0.0;
double A5=0.0;
#endif

#if _SlipTop
scalar slipwall[];
#endif 
 
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
    h         = atof (argv[6]);
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
  #if _MD
  if (argc > 12)
    sign_lev  = atoi (argv[12]);
  #endif  
  #if _SlipTop
  if (argc > 13)
    yh        = atof (argv[13]);
  #endif
  #if _OscilDrop
  if (argc > 14)
    A2        = atof (argv[14]);
  if (argc > 15)
    A3        = atof (argv[15]);
  #endif

  size (L);
  init_grid (256);

  /**
  The dimensionless parameters can then be computed based on 
  rho_l, sigma, and D. */

  rho1 = 1., rho2 = rhog; 
  mu1 = mul, mu2 = mug;
  f.sigma = 1.;
    fprintf(fout,"checkpoint 0 depth=%d\n", max_level); 
 
#if 0 
  gx  = gravity*1100.*sq(2.7e-3)/0.0483;
#endif 
  //TOLERANCE = 1.e-5; 
  //NITERMIN = 5;
  run();
}

#if _MD
/**
Incorporate manifold death
*/
scalar thin_structure[];

event calc_and_print (t = 0.521){
/*  thin_structure.prolongation = myprolongation;
  thin_structure.restriction = myrestrict;
  thin_structure.refine = refine_injection;
*/
  fprintf(stdout,"t=%g\n",t); 
    fprintf(fout,"checkpoint 0 l_sign=%d, depth=%d\n",sign_lev, depth()); 
  foreach(){
    phii[] = 2*f[] - 1;
    sign[] = 7;
   // fprintf(fout,"checkpoint 1\n"); 
  }  
    /** We choose at which level (`l_sign`) we want to compute the quadratic moments and the signature. Then, the indicator function `phii` is restricted at level `l_sign`.*/ 
  int l_sign = sign_lev;
  
    //fprintf(fout,"checkpoint 1\n"); 
  for (int ilev = depth() - 1; ilev >= l_sign; ilev--){ 
  //for (int ilev = max_level - 1; ilev >= l_sign; ilev--){ 
    foreach_level(ilev){
      if(is_refined(cell))
      restriction_average(point, phii);
    }
  }
    printf("checkpoint 2 l_sign=%d\n",l_sign); 
  compute_signature_neigh_level (f, phii, sign, l_sign);
  
  /** 
   The signature `sign` is available only at the level `l_sign`. 
   We need to prolong it onto the finest grid. */
  if (pid()==0){ 
    printf("time %g level used for moments %d and depth is %d \n", t, l_sign, depth()); 
    printf("checkpoint 3\n");}
  for (int ilev = l_sign; ilev < depth(); ilev++){
  //for (int ilev = l_sign; ilev < max_level; ilev++){
    foreach_level(ilev){
      sign.prolongation = phii.prolongation = refine_injection;
      if(is_refined(cell)){
        sign.prolongation (point, sign);
        phii.prolongation (point, phii);
      }
    }
  }
/*
  foreach(){
	  if((f[] > 0.01 && sign[] == -1) ||(f[] < 0.99 && sign[] == 2))
		  sign[] = 1;
  }
  foreach(){
      if(sign[]==1){
	      thin_structure[] = 1;
      }else{
          thin_structure[] = 0;
      }
  }
  boundary({thin_structure});
  restriction({thin_structure});
*/int count=0;
/*  foreach_level(l_sign){
      fprintf(fout, "f=%g, sign=%g, count=%d\n",f[],sign[],count);
      count++;
  }*/
  count =0;
  foreach(){
      if(f[]>0.99 && f[] <1.01){
      //printf("f=%g,count=%d\n",f[],count);
      }
      count++;
  }
  //printf("number of cells=%d\n",count);
  change_topology (f, sign, M, l_sign);
  //change_topology (f, sign, M, l_sign, max_change, large);  
  /**if(ON==1){
      change_topology (f, sign, M, l_sign, &ON);
  }
*/
/** Output the result in HDF5 format*/
 /** scalar * list = {sign,f};
  vector * vlist = {u};
  output_ppm (sign, linear = true, file = "f.mp4", n=200);
*/
}
#endif

/**
The initial drop is spherical. */

event init (t = 0)
{ list_print(all,stdout);
  if (!restore (file = "dump")) {

    double x0 = 0.5+h; 
    double Rd = 0.5; 

    refine (  sq(y)+sq(x-x0)<sq(0.6) && sq(y)+sq(x-x0)>sq(0.4) && level < max_level );
#if _OscilDrop
    fraction (f, - sq(y) - sq(x - x0)
                 + sq(Rd*(1.  + A2* 1./2.*(  3.*pow((x-x0),2.)/pow((sq(x-x0)+sq(y)),1.)
                                          -   1.)
                              + A3* 1./2.*(  5.*pow((x-x0),3.)/pow((sq(x-x0)+sq(y)),1.5)
                                          -   3.*pow((x-x0),1.)/pow((sq(x-x0)+sq(y)),0.5))
                              + A4* 1./8.*( 35.*pow((x-x0),4.)/pow((sq(x-x0)+sq(y)),2.)
                                          -  30.*pow((x-x0),2.)/pow((sq(x-x0)+sq(y)),1.)
                                          +   3.)
                              + A5* 1./8.*( 63.*pow((x-x0),5.)/pow((sq(x-x0)+sq(y)),2.5)
                                          -  70.*pow((x-x0),3.)/pow((sq(x-x0)+sq(y)),1.5)
                                          +  15.*pow((x-x0),1.)/pow((sq(x-x0)+sq(y)),0.5)))));
#else 
    fraction (f, -sq(y)-sq(x-x0)+sq(Rd));
#endif 

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

double xc;
event logfile (i+=10)
{
  scalar posy[],posx[],posx_y0[],udrop[],xdrop[];
  position (f,posx,{1,0});
  position (f,posy,{0,1});
  position (f,posx_y0,{1,0});

  double area=0.,vol=0.,ke=0.,ud=0.,xd=0.,s_area=0;

  foreach(reduction(+:area) reduction(+:vol) reduction(+:ke)
          reduction(+:ud) reduction(+:xd)) {
    #if 1 
    if ( f[] <= 1e-6 || f[] >= 1. - 1e-6 ) 
    {
      posx[] = nodata;
      posy[] = nodata;
      posx_y0[] = nodata;
    }
    #endif 
    if ( y<-Delta || y > Delta ) {
      posx_y0[] = nodata;
    }

    /** statistics in axisymmetric geometry */
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      /** interfacial area */
      coord n = interface_normal (point, f), p;
      double alpha = plane_alpha (f[], n);
      area += pow(Delta, 1.)*plane_area_center (n, alpha, &p)*2.*pi*posy[];
    }
    
    double dv_axi = pow(Delta, 2.)*2.*pi*y;

    /** Volume */
    if (f[] > 1e-6 ) {
      vol += dv_axi*f[];

      /** kinetic energy  */
      foreach_dimension() {
        ke += dv_axi*sq(u.x[]);
      }

      /** mean velocity*/
      ud += dv_axi*f[]*u.x[];

      /** centroid */
      xd += dv_axi*f[]*x;
    }
  }
  ke /= 2.;
  ud /= vol;
  xd /= vol;

  stats sx = statsf (posx);
  stats sy = statsf (posy);
  stats sx_y0 = statsf (posx_y0);

  s_area = interface_area(f);
  
  if ( i == 0 )
    fprintf(ferr,
      "#1: t; 2: dt; 3: xc; 4: uc; 5:y_max; 6:x_max; 7: x_min; 8: x_y0_max; 9: x_y0_min; 10: vol 11: area 12: KE; 13: n_grid; 14: cput; 15: speed\n");

  if ( i >10 ) 
    fprintf (ferr, "%g %g %g %g %g %g %g %g %g %g %g %g %ld %g %g \n",t, dt, xd, ud, sy.max, sx.max, sx.min, sx_y0.max, sx_y0.min, vol, s_area, ke, grid->tn, perf.t, perf.speed);

xc = xd;

}

#if !_MPI
/** output snapshot of simulation */
event interface (t += t_out) {

  char name[80];
  sprintf (name, "infc-%05.3f.dat", t);
  FILE * fp1 = fopen (name, "w");
  output_facets (f,fp1);

}
#endif 
  
event snapshot (t += t_out; t<=t_end ) {
  char name[80];
#if 0
  sprintf (name, "snapshot-%05.3f.gfs", t);
  output_gfs (file = name, t=t, list = {f,u,p});
#endif

  sprintf (name, "dump-%05.3f", t);
  p.nodump = false;
  dump (file = name); // so that we can restart
}


#if _RemoveDrop

double vol_cut = 0.001;// we remove droplets with volume smaller than vol_cut
event remove_drops( i+=10)
{
    scalar m[];
    foreach()
	m[] = f[] > 0;
    int n = tag(m);
    double v[n];
    int remove_flag[n];

    fprintf(stdout,"number of drops=%d\n",n);

    for(int j = 0; j < n; j++){
        v[j] = 0;
	remove_flag[j] = 0;
    }
    foreach_leaf()
	if(m[] > 0){
	    int j = m[]-1;
	    v[j] += dv()*f[];
	}

#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

    for(int j = 0; j < n; j++)
	if(v[j] < vol_cut)
            remove_flag[j] = 1;
    foreach()
	if(m[] > 0){
	    int j = m[] - 1;
	    if(remove_flag[j]==1)
		f[] = 0;
	}
}

#endif



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

#if _UNREFINE
  unrefine ( x>0.5*L && level >= max_level-2 );
#endif
}

/**
We add the acceleration of gravity (unity) in x direction. */

#if 0 
event acceleration (i++) {
  face vector av = a;
  foreach_face(x)
    av.x[] = gx;
}
#endif

#if _SlipTop
event static_slipwall (i++) {
  fraction (slipwall, y-yh);

  foreach()
    u.y[] = (1.-slipwall[])*u.y[];
  boundary (u.y);
}
#endif

#if 1
event runtime (i += 10) {
//  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
/*    if (perf.t/60 >= maxruntime) {
        dump (file = "dump"); // so that we can restart
        return 1; // exit
    }
    */
    if(xc>L-2){//we will stop if drop is 2 diameter away from the right wall
        dump (file = "dump"); // so that we can restart
        return 1; // exit
    }

}
#endif
