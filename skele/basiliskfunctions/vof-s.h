#ifndef  _skel
#define _skel = 1
#endif


int doSmooth = 1;
#if dimension == 2
double *getBounds(double x, double y, double Delta){
  double *retpoint = malloc(4*sizeof(double));//[Xmin,Xmax,Ymin,Ymax]
#else
double *getBounds(double x, double y,double z, double Delta){
  double *retpoint = malloc(6*sizeof(double));
#endif
  retpoint[0] = x-Delta/2;
  retpoint[1] = x+Delta/2;
  retpoint[2] = y-Delta/2;
  retpoint[3] = y+Delta/2;
#if dimension == 3
  retpoint[4] = z-Delta/2;
  retpoint[5] = z+Delta/2;
#endif
  return retpoint;
}

int pointInsideCell(double *bounds, double *point){
  if(bounds[0] < point[0] && bounds[1] > point[0]){
    if(bounds[2] < point[1] && bounds[3] > point[2]){
#if dimension == 3
      if(bounds[4] < point[2] && bounds[5] > point[2])
#endif
        return 1;
    }
  }
  return 0;
}

#if _MPI
int comm_size = 0;    
//useful MPI related comunication for skeleton/vof related data :)
//NOTE:
//we define our list of points as a double**; this will allow us to alloc a total size of the entire field
//Only at the top level, ie points = malloc(totallength*sizeof(double*))
//Then when loading in a point we can alloc the full space we will need for said point ie [x,y,z,r]
//This will help with keeping consistent id's and refrences through mpi; and allow for easy sorting
//and transfer of points which will be relevant in each rank
int *gatherLengthMPI(int inlen){
  int *sumarr = malloc(comm_size*sizeof(int));
  sumarr[pid()] = inlen;
  for(int i = 0; i < comm_size; i++){
    if(i == pid()){
      //i will recieve from all pids
      for(int j = 0; j < comm_size; j++){
        if(j != i){
          MPI_Recv(&sumarr[j],1,MPI_INT,j,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
      }
    }
    else{
      MPI_Send(&inlen,1,MPI_INT,i,0,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return sumarr;
}
#endif 

//structure assoicated with 'fid', holds skeletons and other useful data in the fields which cant be held in scalars
struct skeledata{
  //relevant position and ID ranges
  double **sPts;
  int *lstart;//local start
  int *lend;//local end
  int *length;//total length
  
  //ID positions neighbors and their respective amounts
  int **neighbors;//position of connected skeletons
  int *nlengths;//position of connected skeletons

  //for holding points if structure exitsts,these get cleaned out each time 
  double **holdpoints;
  int *hplength;
};
void addneighbor(struct skeledata *skel,int nref0,int nref1){
  //this function takes in a structure and adds to the neighbor field of it
  printf("neighcalc %d - %d",nref0,nref1);
  skel->neighbors[nref0] = realloc(skel->neighbors[nref0],(skel->nlengths[nref0] + 1)*sizeof(int));
  printf("r0 good\n");
  skel->neighbors[nref0][skel->nlengths[nref0]] = nref1;
  printf("r0set good\n");
  printf("r0l b4 %d\n",skel->nlengths[nref0]);
  skel->nlengths[nref0]++;
  printf("r0l good %d\n",skel->nlengths[nref0]);
  skel->neighbors[nref1] = realloc(skel->neighbors[nref1],(skel->nlengths[nref1] + 1)*sizeof(int));
  printf("r1 good\n");
  skel->neighbors[nref1][skel->nlengths[nref1]] = nref0;
  printf("r1set good\n");
  printf("r1l b4 %d\n",skel->nlengths[nref1]);
  skel->nlengths[nref1]++;
  printf("r1l good %d\n",skel->nlengths[nref1]);
}

//We use a structure associated with transportfid, this allows us to have multiple normals and interfaces within a cell, and still be used
struct reconstructTransport{
  int **interfaceLocations;//holds the ids of interfaces, eg[[0,1],[5,8],[3],[6,2],..], uses skeledata id as location 
  int *interfacecount;//holds amount of interfaces which exist at the curent location, eg[2,2,1,2,...]
  int *illength;//length of finterfaceLocatonis and count
  double **normals;//individual normals
  double *alphas;//individual alphas
  int *length;//length of normals and alphas
};
void deleteRT(struct reconstructTransport *RTD){
  if(RTD != NULL){
    for(int i = 0; i < *RTD->length; i++){
      free(RTD->normals[i]);
    }
    free(RTD->alphas);
    free(RTD->normals);
    free(RTD->length);
    for(int i = 0; i < *RTD->illength; i++){
        free(RTD->interfaceLocations[i]);
    }
    free(RTD->interfaceLocations);
    free(RTD->interfacecount);
    free(RTD->illength);
    free(RTD);
  }
}

attribute {
  scalar * tracers, c, fid,rid;
  bool inverse;
  struct skeledata *skel;
  struct reconstructTransport *recont;
}

#include "fractions.h"
#include "skele/skeleton.h"

extern scalar * interfaces;
extern scalar * interfacefid;

extern face vector uf;
extern double dt;

extern double max_level;//max level of sim
extern double L;//domain size

//for ease of computations on MPI we track active PID's with interfaces, and block out rest of processors
#if _MPI
int *active_PID;
#endif

/**
The gradient of a VOF-concentration `t` is computed using a standard
three-point scheme if we are far enough from the interface (as
controlled by *cmin*), otherwise a two-point scheme biased away from
the interface is used. */

foreach_dimension()
static double vof_concentration_gradient_x (Point point, scalar c, scalar t)
{
  static const double cmin = 0.5;
  double cl = c[-1], cc = c[], cr = c[1];
  if (t.inverse)
    cl = 1. - cl, cc = 1. - cc, cr = 1. - cr;
  if (cc >= cmin && t.gradient != zero) {
    if (cr >= cmin) {
      if (cl >= cmin) {
	if (t.gradient)
	  return t.gradient (t[-1]/cl, t[]/cc, t[1]/cr)/Delta;
	else
	  return (t[1]/cr - t[-1]/cl)/(2.*Delta);
      }
      else
	return (t[1]/cr - t[]/cc)/Delta;
    }
    else if (cl >= cmin)
      return (t[]/cc - t[-1]/cl)/Delta;
  }
  return 0.;
}

/**
On trees, VOF concentrations need to be refined properly i.e. using
volume-fraction-weighted linear interpolation of the concentration. */

#if TREE
static void vof_concentration_refine (Point point, scalar s)
{
  scalar f = s.c;
  if (cm[] == 0. || (!s.inverse && f[] <= 0.) || (s.inverse && f[] >= 1.))
    foreach_child()
      s[] = 0.;
  else {
    coord g;
    foreach_dimension()
      g.x = Delta*vof_concentration_gradient_x (point, f, s);
    double sc = s.inverse ? s[]/(1. - f[]) : s[]/f[], cmc = 4.*cm[];
    foreach_child() {
      s[] = sc;
      foreach_dimension()
	s[] += child.x*g.x*cm[-child.x]/cmc;
      s[] *= s.inverse ? 1. - f[] : f[];
    }
  }
}
//Because we will need to transport id's around effectively, we define a prolongation and refinement method



/**
On trees, we need to setup the appropriate prolongation and
refinement functions for the volume fraction fields. */
static void fid_prolongation(Point point, scalar fid){
    double val = fid[];
    if(fid[] != 0.){
        //all parent data applies to children data if it doesnt exist
        foreach_child(){
            fid[] = val;
        }
    }
    else{
        //if no parent data, no child data, though this shouldnt happen
        foreach_child(){
           fid[] = 0.;
        }
    }
}

static inline void fid_restriction(Point point, scalar fid){
    double val = 0.;
    int pass = 0.;
    foreach_child(){//we pick the first child to hold a point when we upscale, allows for even distribution of points/smoothening
        if(!pass && fid[] != 0){
            val = fid[];
            pass++;
        }
    }
    fid[] = val;
}

event defaults (i = 0)
{
  //ensure that we give each interface access to all skeledata,vof data, and id data 
  scalar c,fid;
  for (c,fid in interfaces,interfacefid) {
    fid.coarsen = fid.restriction = fid_restriction;
    fid.refine = fid.prolongation = fid_prolongation;
    //fid.dirty = true;
    struct skeledata *sd = malloc(sizeof(struct skeledata));
    sd->length = malloc(sizeof(int));
    sd->lstart = malloc(sizeof(int));
    sd->lend = malloc(sizeof(int));
    *sd->length = 0,*sd->lstart = 0,*sd->lend = 0; 
    c.refine = c.prolongation = fraction_refine;
    c.dirty = true;
    c.fid = fid;
    fid.skel = sd;
    fid.c = c;
    c.skel = sd;
    scalar * tracers = c.tracers;
    for (scalar t in tracers) {
      printf("made inside tracer\n");
      t.restriction = restriction_volume_average;
      t.refine = t.prolongation = vof_concentration_refine;
      t.dirty = true;
      t.c = c;
      t.fid = fid;
      t.skel = sd;
    }
  }
}
#endif // TREE

/**
Boundary conditions for VOF-advected tracers usually depend on
boundary conditions for the VOF field. */

event defaults (i = 0)
{
  for (scalar c in interfaces) {
    scalar * tracers = c.tracers;
    for (scalar t in tracers)
      t.depends = list_add (t.depends, c);
  }
  //Last thing we do is allocate Active_PID if needed
#if _MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  active_PID = calloc(comm_size , sizeof(int));//ensure we have a marker for all active cells 
#endif
}



/**
We need to make sure that the CFL is smaller than 0.5 to ensure
stability of the VOF-skeleton scheme. */


event stability (i++) {
  if (CFL > 0.5)
    CFL = 0.5;
}


/**
## One-dimensional advection

The simplest way to implement a multi-dimensional VOF advection scheme
is to use dimension-splitting i.e. advect the field along each
dimension successively using a one-dimensional scheme.

We implement the one-dimensional scheme along the x-dimension and use
the [foreach_dimension()](/Basilisk C#foreach_dimension) operator to
automatically derive the corresponding functions along the other
dimensions. */

foreach_dimension()
static void sweep_x (scalar c, scalar cc, scalar * tcl)
{
  vector n[];
  scalar alpha[], flux[];
  double cfl = 0.;
 
  /**
  If we are also transporting tracers associated with $c$, we need to
  compute their gradient i.e. $\partial_xf_j = \partial_x(t_j/c)$ or
  $\partial_xf_j = \partial_x(t_j/(1 - c))$ (for higher-order
  upwinding) and we need to store the computed fluxes. We first
  allocate the corresponding lists. */

  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }

    /**
    The gradient is computed using the "interface-biased" scheme above. */

    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	    gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }
  
  /**
  We reconstruct the interface normal $\mathbf{n}$ and the intercept
  $\alpha$ for each cell. Then we go through each (vertical) face of
  the grid. */

  reconstruction (c, n, alpha);
  foreach_face(x, reduction (max:cfl)) {

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.
    
    ![Volume fraction flux](figures/flux.svg)
    
    When the upwind cell is entirely full or empty we can avoid this
    computation. */

    double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
      rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
			  (coord){-0.5, -0.5, -0.5},
			  (coord){s*un - 0.5, 0.5, 0.5});
    
    /**
    Once we have the upwind volume fraction *cf*, the volume fraction
    flux through the face is simply: */

    flux[] = cf*uf.x[];

    /**
    If we are transporting tracers, we compute their flux using the
    upwind volume fraction *cf* and a tracer value upwinded using the
    Bell--Collela--Glaz scheme and the gradient computed above. */
    
    scalar t, gf, tflux;
    for (t,gf,tflux in tracers,gfl,tfluxl) {
      double cf1 = cf, ci = c[i];
      if (t.inverse)
	cf1 = 1. - cf1, ci = 1. - ci;
      if (ci > 1e-10) {
	double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
	tflux[] = ff*cf1*uf.x[];
      }
      else
	tflux[] = 0.;
    }
  }
  delete (gfl); free (gfl);
  
  /**
  We warn the user if the CFL condition has been violated. */

  if (cfl > 0.5 + 1e-6)
    fprintf (ferr, 
	     "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
	     cfl - 0.5), fflush (ferr);

  /**
  Once we have computed the fluxes on all faces, we can update the
  volume fraction field according to the one-dimensional advection
  equation
  $$
  \partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  $$
  The first term is computed using the fluxes. The second term -- which is
  non-zero for the one-dimensional velocity field -- is approximated using
  a centered volume fraction field `cc` which will be defined below. 

  For tracers, the one-dimensional update is simply
  $$
  \partial_tt_j = -\nabla_x\cdot(\mathbf{u}_f t_j)
  $$
  */

#if !EMBED
  foreach() {
    c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
    scalar t, tc, tflux;
    for (t, tc, tflux in tracers, tcl, tfluxl)
      t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  }
#else // EMBED
  /**
  When dealing with embedded boundaries, we simply ignore the fraction
  occupied by the solid. This is a simple approximation which has the
  advantage of ensuring boundedness of the volume fraction and
  conservation of the total tracer mass (if it is computed also
  ignoring the volume occupied by the solid in partial cells). */
  
  foreach()
    if (cs[] > 0.) {
      c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
      scalar t, tc, tflux;
      for (t, tc, tflux in tracers, tcl, tfluxl)
	t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
    }
#endif // EMBED

  delete (tfluxl); free (tfluxl);
}

/**
## Multi-dimensional advection

The multi-dimensional advection is performed by the event below. */


void skelPushValues(double *val1, double *val2,int amount){
    //for setting all values of 1 equal to value 2
    for(int i = 0; i < amount; i++){
        val1[i] = val2[i];
    }
}
#include "fractions-s.h"

//Here wee have the 1D advection implementation for vof-skeleton combination
foreach_dimension()
static void skeladv_x (scalar c, scalar cc, scalar * tcl)
{
  vector n[];
  scalar rid[];//rid is the flag for if an interface has multiple within a cell 
  scalar alpha[], flux[];
  alpha.rid = rid;
  //n.x.rid = rid;
  //scalar fid = c.fid;
  //struct skeledata *sd = c.skel;
  //set up our structure needed for appropriate calculation
  struct reconstructTransport *rt = NULL;
  alpha.recont = rt;
  //n.x.recont = rt;
  //build reconstruction id, based on if skeleton points are present near, for all cells
  double cfl = 0.;
  scalar * tracers = c.tracers, * gfl = NULL, * tfluxl = NULL;
  if (tracers) {
    for (scalar t in tracers) {
      scalar gf = new scalar, flux = new scalar;
      gfl = list_append (gfl, gf);
      tfluxl = list_append (tfluxl, flux);
    }
    foreach() {
      scalar t, gf;
      for (t,gf in tracers,gfl)
	    gf[] = vof_concentration_gradient_x (point, c, t);
    }
  }
  //we reconstruction vofcells independantly of our skeleton cells
  reconstruction_v (c, n, alpha);
  reconstruction_s (c, n, alpha);
  //reconstructionCorrectionSkel(c, n, alpha);
  foreach_face(x, reduction (max:cfl)) {

    /**
    To compute the volume fraction flux, we check the sign of the velocity
    component normal to the face and compute the index `i` of the
    corresponding *upwind* cell (either 0 or -1). */

    double un = uf.x[]*dt/(Delta*fm.x[] + SEPS), s = sign(un);
    //int i = -(s + 1.)/2.;

    /**
    We also check that we are not violating the CFL condition. */

#if EMBED
    if (cs[] >= 1.)
#endif
    if (un*fm.x[]*s/(cm[] + SEPS) > cfl)
      cfl = un*fm.x[]*s/(cm[] + SEPS);

    /**
    If we assume that `un` is negative i.e. `s` is -1 and `i` is 0, the
    volume fraction flux through the face of the cell is given by the dark
    area in the figure below. The corresponding volume fraction can be
    computed using the `rectangle_fraction()` function.
    
    ![Volume fraction flux](figures/flux.svg)
    
    When the upwind cell is entirely full or empty we can avoid this
    computation. */

  //  double cf = (c[i] <= 0. || c[i] >= 1.) ? c[i] :
  //    rectangle_fraction ((coord){-s*n.x[i], n.y[i], n.z[i]}, alpha[i],
  //  		  (coord){-0.5, -0.5, -0.5},
  //  		  (coord){s*un - 0.5, 0.5, 0.5});
  //  
  //  /**
  //  Once we have the upwind volume fraction *cf*, the volume fraction
  //  flux through the face is simply: */

  //  flux[] = cf*uf.x[];

  //  /**
  //  If we are transporting tracers, we compute their flux using the
  //  upwind volume fraction *cf* and a tracer value upwinded using the
  //  Bell--Collela--Glaz scheme and the gradient computed above. */
  //  
  //  scalar t, gf, tflux;
  //  for (t,gf,tflux in tracers,gfl,tfluxl) {
  //    double cf1 = cf, ci = c[i];
  //    if (t.inverse)
  //  cf1 = 1. - cf1, ci = 1. - ci;
  //    if (ci > 1e-10) {
  //  double ff = t[i]/ci + s*min(1., 1. - s*un)*gf[i]*Delta/2.;
  //  tflux[] = ff*cf1*uf.x[];
  //    }
  //    else
  //  tflux[] = 0.;
  //  }
  }
  delete (gfl); free (gfl);
  //
  ///**
  //We warn the user if the CFL condition has been violated. */

  //if (cfl > 0.5 + 1e-6)
  //  fprintf (ferr, 
  //       "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n", 
  //       cfl - 0.5), fflush (ferr);

  ///**
  //Once we have computed the fluxes on all faces, we can update the
  //volume fraction field according to the one-dimensional advection
  //equation
  //$$
  //\partial_tc = -\nabla_x\cdot(\mathbf{u}_f c) + c\nabla_x\cdot\mathbf{u}_f
  //$$
  //The first term is computed using the fluxes. The second term -- which is
  //non-zero for the one-dimensional velocity field -- is approximated using
  //a centered volume fraction field `cc` which will be defined below. 

  //For tracers, the one-dimensional update is simply
  //$$
  //\partial_tt_j = -\nabla_x\cdot(\mathbf{u}_f t_j)
  //$$
  //*/

#if !EMBED
  //foreach() {
  //  c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  //  scalar t, tc, tflux;
  //  for (t, tc, tflux in tracers, tcl, tfluxl)
  //    t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/(cm[]*Delta);
  //}
#else // EMBED
  ///**
  //When dealing with embedded boundaries, we simply ignore the fraction
  //occupied by the solid. This is a simple approximation which has the
  //advantage of ensuring boundedness of the volume fraction and
  //conservation of the total tracer mass (if it is computed also
  //ignoring the volume occupied by the solid in partial cells). */
  //
  //foreach()
  //  if (cs[] > 0.) {
  //    c[] += dt*(flux[] - flux[1] + cc[]*(uf.x[1] - uf.x[]))/Delta;
  //    scalar t, tc, tflux;
  //    for (t, tc, tflux in tracers, tcl, tfluxl)
  //      t[] += dt*(tflux[] - tflux[1] + tc[]*(uf.x[1] - uf.x[]))/Delta;
  //  }
#endif // EMBED

  delete (tfluxl); free (tfluxl);
  delete({rid});
  rt = alpha.recont;
  deleteRT(rt);
}


//hold passable dimension for skeleton, is static
void inject_skele(scalar c, double **inputSkeleton, int *length, int reali){
  //This function will input ours list of skeletons [x,y,r] or [x,y,z,r]  which are [p0..pl] long into appropiate scalars
  //load out id field and skeleton field
  scalar fid = c.fid;
  struct skeledata *skel = c.skel;
  int *addid = calloc(*length,sizeof(int));
  double **computeSkel = malloc(*length * sizeof(double*));
  int csindx = 0;
  for(int i = 0; i < *length; i++){
    Point point;
    if(!addid[i]){//only add if not calculated yet
      addid[i]++;
#if dimension == 2
      point = locate(inputSkeleton[i][0],inputSkeleton[i][1]);
      double *bounds = getBounds(x,y,Delta);
      if(fid[] < 0.99){//ensure uncalculated region
        computeSkel[csindx] = calloc(3,sizeof(double));//alloc for a new point
#else
      point = locate(inputSkeleton[i][0],inputSkeleton[i][1],inputSkeleton[i][2]);
      double *bounds = getBounds(x,y,z,Delta);
      if(fid[] < 0.99){//ensure uncalculated region
        computeSkel[csindx] = calloc(4,sizeof(double));//alloc for a new point
#endif
        int q, tot =  0;
        for(q = 0; q < dimension+1; q++){
          computeSkel[csindx][q] = inputSkeleton[i][q];
        }
        tot++;
        for(int j = i+1; j < *length;j++){
          //secondary loop for averaging
          if(!addid[j] && pointInsideCell(bounds,inputSkeleton[j])){
            for(q = 0; q < dimension+1; q++){
              computeSkel[csindx][q] += inputSkeleton[j][q];
            }
            tot++;
            addid[j]++;
          }
        }
        for(q = 0; q < dimension+1; q++){
          computeSkel[csindx][q] = computeSkel[csindx][q] / tot;
        }
        fid[] = (double)csindx + 1.;//we build ID's up and set local counts 
        csindx++;
      }
      else{
          printf("error building location id already exists! %f\n",fid[]);
      }
      free(bounds);
    }
  }
  //now that we are done with our input skeleton, we free it up
  if(inputSkeleton != NULL){
    for(int q = 0; q < *length; q++){
      free(inputSkeleton[q]);
    }
    free(inputSkeleton);
    inputSkeleton = NULL;
  }
#if _MPI
  int istart = 0;
  int iend = 0;
  int ilength = 0;
  int *mpilengths = gatherLengthMPI(csindx);
  for(int i = 0; i < comm_size; i++){
    if(i <= pid()){
      iend += mpilengths[i];//makes our ending id at end of local data 
    }
    if(i < pid()){
      istart += mpilengths[i];//start adjusts to our real start
    }
    ilength += mpilengths[i];//all have same length
  }
  //printf("%d - start(%d) end(%d) total(%d)\n",pid(),istart,iend,ilength);
  //for(int i = 0; i < csindx; i++){
  //  printf("%d has found point%d :[%f,%f,%f]\n",pid(),i,computeSkel[i][0],computeSkel[i][1],computeSkel[i][2]);
  //}
  free(mpilengths);
#else
  int istart = 0;
  int iend=csindx;
  int ilength=csindx;
#endif
  if(ilength > 0){
    //now we can set our struct values and assign id's
    //NOTE: our positions will be shifted up one; this is beacuse id = 0 means 
    *skel->length = ilength;
    *skel->lstart = istart;
    *skel->lend = iend;
    //Firstly we want to provide ample space to gather 
    //points and keep consisten numbering;
    double **pointfull = malloc(ilength * sizeof(double*));//holds our points in relevant positions
    int **neighborfull = malloc(ilength * sizeof(double*));//holds our neighbor points to relevant positions
    int *neighborfulllength = calloc(ilength , sizeof(int));//holds the amount of neighbors which we have
    for(int i = 0; i < ilength; i++){
        pointfull[i] = calloc(dimension + 1,sizeof(double));//set space for dim + r
        neighborfull[i] = calloc(1,sizeof(double));//only set space for one for now...
    }
#if _MPI
    int *inow = malloc(comm_size*sizeof(int));
    inow[pid()] = istart;
    foreach(reduction(+:inow[:comm_size])){//first we gather everything within our current region
      if(fid[] > 0.99){
        int localtag = (int)fid[] - 1;//currently holds our position of computeSkel
        //Now we want to transfer our allocated computeSkel into our structure
        skelPushValues(pointfull[inow[pid()]],computeSkel[localtag],dimension + 1);//swap compute to our pointfull
        fid[] = (double)inow[pid()] + 1.;//apply upshift
        inow[pid()]++;
      }
    }
    free(inow);//clean up
#else
    int inow = istart;
    foreach(){//first we gather everything within our current region
      if(fid[] != 0.){
        int localtag = (int)fid[] - 1;//currently holds our position of computeSkel
        //Now we want to transfer our allocated computeSkel into our structure
        skelPushValues(pointfull[inow],computeSkel[localtag],dimension + 1);//swap compute to our pointfull
        fid[] = (double)inow + 1.;//apply upshift
        inow++;
      }
    }
#endif
    skel->sPts = pointfull;
    skel->neighbors = neighborfull;
    skel->nlengths = neighborfulllength;
    //now that we have put our values safely into their own space we get rid of other allocations
    for(int i = 0; i < *length; i++){
      if(computeSkel[i] != NULL)free(computeSkel[i]);
    }
    //Next we can go through and add in our intermediate values
    boundary({fid});
    for(int q = max_level - 1; q >= 0; q--){
        //printf("%drunning bound @ %d; treefull? %d\n",pid(),q,tree_is_full());
        boundary_level({fid},q);
    }
    foreach(){
      if(fid[] == 0.){
        double val = 0.;
        foreach_neighbor(1){
          if(fid[] > 0.99){
            //printf("val1 = %f -> %f\n",val,fid[]);
            val = -1. * fid[];
            //printf("val2 = %f -> %f\n",val,fid[]);
          }
        }
        fid[] = val;
      }
      if(fid[] > 0.99){
        int nref = (int)fid[] - 1;//shit down
        double href = fid[];
        foreach_neighbor(1){
          //if neighbor has val we add it to our structure
          if(fid[] > 0.99 && fid[] != href){
            printf("adding %d & %f\n",nref,fid[]);
            addneighbor(skel,nref,(int)fid[] - 1);
          }
        }
      }
    }
  }
  free(computeSkel);
  free(addid);
  //at the end we want to be sure and clean each
}

double salpha = 40 * PI / 180;/////INPUT ANGLE TO SKELETON 
//scan skel will detect & create skeletons
void scanSkel(scalar * interfaces,int i){
  for (scalar c in interfaces) {
    scalar fid = c.fid;
    boundary({fid});
    for(int q = max_level - 1; q >= 0; q--){
      //printf("%drunning bound @ %d; treefull? %d\n",pid(),q,tree_is_full());
      boundary_level({fid},q);
    }
    
    struct skeledata *skel = c.skel;
    if(*skel->length == 0){//what to do when we have no structure defined yet
#if _MPI
      foreach(){
        if(!active_PID[cell.pid]){//force a faster blockoff when itterating bc basilisk no break rule
          if(c[] > 1e-6 && c[] < 1.-1e-6){
            active_PID[cell.pid] = 1;
          }
        }
      }
#endif
      double **calcskeleton = NULL;//skeleton container
      int skelelength = 0;
      printf("calc @ %d\n",i);
#if _MPI
      calcSkeletonMPI(c,&salpha,max_level,L,(double)i/1000.,&calcskeleton,&skelelength,doSmooth);
#else
      calcSkeleton(c,&salpha,max_level,L,(double)i/1000.,&calcskeleton,&skelelength,doSmooth);
#endif
      printf("sl=%d\n",skelelength);
      //inject_skele(c,calcskeleton,&skelelength,i);//insert skeletons into field
    }
  }

}

void real_vof_advection (scalar * interfaces, int i)
{ 
  for (scalar c in interfaces) {
    /**
    We first define the volume fraction field used to compute the
    divergent term in the one-dimensional advection equation above. We
    follow [Weymouth & Yue, 2010](/src/references.bib#weymouth2010) and use a
    step function which guarantees exact mass conservation for the
    multi-dimensional advection scheme (provided the advection velocity
    field is exactly non-divergent). */
    //scalar fid = c.fid;
    struct skeledata *sd = c.skel;
    scalar cc[], * tcl = NULL, * tracers = c.tracers;    
    for (scalar t in tracers) {
      scalar tc = new scalar;
      tcl = list_append (tcl, tc);
#if TREE
      if (t.refine != vof_concentration_refine) {
	    t.refine = t.prolongation = vof_concentration_refine;
	    t.restriction = restriction_volume_average;
	    t.dirty = true;
	    t.c = c;
      }
#endif // TREE
    }
    foreach() {
      cc[] = (c[] > 0.5);
      scalar t, tc;
      for (t, tc in tracers, tcl) {
	    if (t.inverse)
	      tc[] = c[] < 0.5 ? t[]/(1. - c[]) : 0.;
	    else
	      tc[] = c[] > 0.5 ? t[]/c[] : 0.;
      }
    }

    /**
    We then apply the one-dimensional advection scheme along each
    dimension. To try to minimise phase errors, we alternate dimensions
    according to the parity of the iteration index `i`. */
    //we determine if we do vof advection or vof-skeleton advection based on if skeletons exist
    if(*sd->length){
      void (* sweep[dimension]) (scalar, scalar, scalar *);
      int d = 0;
      foreach_dimension()
        sweep[d++] = skeladv_x;
      for (d = 0; d < dimension; d++)
        sweep[(i + d) % dimension] (c, cc, tcl);
      delete (tcl), free (tcl);
    }
    else{
      void (* sweep[dimension]) (scalar, scalar, scalar *);
      int d = 0;
      foreach_dimension()
        sweep[d++] = sweep_x;
      for (d = 0; d < dimension; d++)
        sweep[(i + d) % dimension] (c, cc, tcl);
      delete (tcl), free (tcl);
    }
  }
}

void vof_advection(scalar * interfaces, int i){
  //new vof advection function will choose what steps to take
  scanSkel(interfaces,i);
  real_vof_advection(interfaces,i);
}

void destroySkel(scalar *interfaces){
  scalar c[];
  for (scalar c in interfaces){
    struct skeledata *skel = c.skel; 
    if(*skel->length > 0){
      for(int i = 0; i < *skel->length; i++){
        //frees any needed lists
        free(skel->sPts[i]);
        free(skel->neighbors[i]);
      }
      free(skel->sPts);
      free(skel->neighbors);
      free(skel->nlengths);
    }
    free(skel->length);
    free(skel->lstart);
    free(skel->lend);
    free(skel);
    skel = NULL;
  }
#if _MPI
  free(active_PID);
#endif
}

event vof (i++){
  vof_advection (interfaces, i);
}
event finalclean (t = end){
  destroySkel(interfaces);
}






/**
## References

~~~bib
@Article{lopez2015,
  title = {A VOF numerical study on the electrokinetic effects in the 
           breakup of electrified jets},
  author = {J. M. Lopez-Herrera and A. M. Ganan-Calvo and S. Popinet and
            M. A. Herrada},
  journal = {International Journal of Multiphase Flows},
  pages = {14-22},
  volume = {71},
  year = {2015},
  doi = {doi.org/10.1016/j.ijmultiphaseflow.2014.12.005},
  url = {http://gfs.sf.net/papers/lopez2015.pdf}
}
~~~
*/
