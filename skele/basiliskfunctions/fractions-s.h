#include "fractions.h"
#include "vof-s.h"
#include "skele/skeleton.h"


struct reconstructIntermediate{
  double **normal;//[nx[0],ny[0],...]
  int **npairs;//holds nref pairs id's incase needed
  int *nlength;
  int **decoder;
  int *dlength;//holds how many refrences each has
#if _MPI
  int *avoidMPI;//keep local value in a foreach...
#endif
  //section for carrying normals out
  double **ncarry;
  double **npoint;
  int *ncarryl;
};


//methods for averaging norms
void splitNorms(double **normsin,int nl,double **norms0,int *n0l, double **norms1,int *n1l,int clean){
  //divide a list of normals into 2 lists
  double *rootn = normsin[0];
  for(int i = 0; i < nl; i++){
    printf("norm in %d => [%f,%f]\n",i,normsin[i][0],normsin[i][1]);
  }
  norms0 = malloc(sizeof(double*));
  norms0[0] = rootn;
  *n0l++;
  for(int i = 1; i < nl; i++){
    //sort normals based on root
    double dp = 0.;
    for(int j = 0; j < dimension; j++){
      dp = dp + rootn[j] * normsin[i][j];//calc dot product
    }
    if(dp > 0.){
      norms0 = realloc((*n0l + 1) * sizeof(double*));
      norms0[*n0l] = normsin[i];
      *n0l++;
    }
    else{
      norms1 = realloc((*n1l + 1) * sizeof(double*));
      norms1[*n1l] = normsin[i];
      *n1l++;
    }
  }
  if(clean){
    //if clean we clean normsin 
    free(normsin);//we dont release points as they are in norms0 and norms1
  }
  for(int i = 0; i < *n0l; i++){
    printf("norm0 out %d => [%f,%f]\n",i,norms0[i][0],norms0[i][1]);
  }
  for(int i = 0; i < *n1l; i++){
    printf("norm1 out %d => [%f,%f]\n",i,norms1[i][0],norms1[i][1]);
  }
}
void averageNormList(double **norm0,int n0l,double **norm1,int n1l,double *nout,int clean){
  //This function will aveage each norm list & return their average norms
  nout = malloc(2*dimension*sizeof(double));
  for(int i = 0; i < n0l; i++){
    for(int j = 0; j < dimension; j++){
      nout[j] = nout[j] + norm0[i][j];
    }
  }
  for(int j = 0; j < dimension; j++){
    nout[j] = nout[j] / n0l;
  }
  for(int i = 0; i < n1l; i++){
    for(int j = 0; j < dimension; j++){
      nout[j+dimension] = nout[j+dimension] + norm1[i][j];
    }
  }
  for(int j = 0; j < dimension; j++){
    nout[j+dimension] = nout[j+dimension] / n1l;
  }
  if(clean){
    for(int i = 0; i < n0l; i++){
      free(norm0[i]);
    }
    for(int i = 0; i < n1l; i++){
      free(norm1[i]);
    }
    free(norm0);
    free(norm1);
  }
}
void gatherSkelNorm(double **innorms, int *inlen,struct reconstructIntermediate *ri,int nref){
  int length = ri->dlength[nref];
  for(int i = 0; i < length; i++){
    innorms = realloc((*inlen + 2)*sizeof(double*));//set room for both normals
    for(int j = 0; j < dimension; j++){
        innorms[*inlen][j]   = ri->normals[ri->decoder[nref][i]][j];
        innorms[*inlen+1][j] = ri->normals[ri->decoder[nref][i]][j+dimension];
    }
    *inlen = *inlen + 2;
  }
}
double *getNormPoint(double *normin,int s,double *skelepoint){
  int start = s ? 0:dimension,end = s ? dimension: 2*dimension;
  double *retnormp = malloc((dimension)*sizeof(double));
  int indx = 0;
  while(start < end){
    retnormp[indx] = skelepoint[indx] + normin[start] * skelepoint[dimension];
    indx++;
    start++;
  }
  return retnormp;
}
//averaging
void normNN(Point point,double *normin,vector n,scalar c){
  //if no neighbors, we need to compute from surounding information
  //We can safely assume that the interface is small here
  //so we gather norms and split them
  double **inputnorms = malloc(10*sizeof(double*)); 
  int inc = 0;
  foreach_neighbor(1){
    if(c[] > 1e-6 && c[] < 1-1e-6){
      inputnorms[inc] = malloc(dimension * sizeof(double));
      int q = 0;
      foreach_dimension()
        inputnorms[inc][q++] = n.x;
      inc++;
    }
  }
  double **norm0,**norm1;
  int n0l=0,n1l=0;
  splitnorms(inputnorms,inc,norm0,&n0l,norm1,&n1l,1);//this will give us norm0&1 and their lengths, first large malloc will be cleaned
  double *nout;
  averageNormList(norm0,n0l,norm1,n1l,nout,1);
  normin = nout;
}
void normON(Point point,double *normin,vector n,scalar c,struct reconstructIntermediate *ri,int nref){
  //if we have one neighbor we will first grab our vofneighbors
  double **inputnorms = malloc(9*sizeof(double*)); 
  int inc = 0;
  foreach_neighbor(1){
    if(c[] > 1e-6 && c[] < 1-1e-6){
      inputnorms[inc] = malloc(dimension * sizeof(double));
      int q = 0;
      foreach_dimension()
        inputnorms[inc][q++] = n.x[];
      inc++;
    }
  }
  gatherSkelNorm(inputnorms,&inc,ri,nref);
  double **norm0,**norm1;
  int n0l=0,n1l=0;
  splitnorms(inputnorms,inc,norm0,&n0l,norm1,&n1l,1);//this will give us norm0&1 and their lengths, first large malloc will be cleaned
  double *nout;
  averageNormList(norm0,n0l,norm1,n1l,nout,1);
  normin = nout;
}
void normALL(Point point,double *normin,scalar c,struct reconstructIntermediate *ri,int nref){
  //if we have one neighbor we will first grab our vofneighbors
  struct skeledata *sd = c.skel;
  double **inputnorms = malloc((2 * sd->nlengths[nref])*sizeof(double*)); 
  int inc = 0;
  gatherSkelNorm(inputnorms,&inc,ri,nref);
  double **norm0,**norm1;
  int n0l=0,n1l=0;
  splitnorms(inputnorms,inc,norm0,&n0l,norm1,&n1l,1);//this will give us norm0&1 and their lengths, first large malloc will be cleaned
  double *nout;
  averageNormList(norm0,n0l,norm1,n1l,nout,1);
  normin = nout;
}
//fractions-s.h holds the needed functions to reconstruct a interface that contains a skeleton
void reconstruction_v(scalar c, vector n, scalar alpha){
  scalar fid = c.fid;
  //similar to default reconstruction, however we only compute on clean vof cells with no skeletons
  foreach(){
    //compute normal for each cell that doesnt have a skeleton point
    if(fid[] < 0.99){
      if (c[] <= 0. || c[] >= 1.) {
        alpha[] = 0.;
        foreach_dimension()
	      n.x[] = 0.;
      }
      else {
        coord m = interface_normal (point, c);
        foreach_dimension()
	      n.x[] = m.x;
        alpha[] = plane_alpha (c[], m);
      }
    }
  }
}
//where reconstruct transport will hold important information about the normals 
//reconstructIntermediate will allow us to build pairs of skeleton points indidually, and then average needed duplicates
void destroyRInter(struct reconstructIntermediate *ri,struct skeledata *sd){
  //responible for cleaning up structure
  for(int i = 0; i < *sd->length; i++){
    free(ri->decoder[i]);
  }
  free(ri->dlength);
  for(int i = 0; i < *ri->nlength; i++){
    free(ri->normal[i]);
    free(ri->npairs[i]);
  }
  free(ri->normal);
  free(ri->npairs);
  free(ri->nlength);
}
//each pair of skeltons should produce 2 normals
//out normals can be bout thorugh trig
//we find the two normals from the two common tangents of the cirlces
//which dont cross their connecting line
void computeSkelNorm(double *skel0,double *skel1,double *output){
  //compute needed values
  int dim = dimension;
  double dist   = getDistance(skel0,skel1,&dim); 
  double dydx   = (skel1[0] - skel0[0]) / (skel1[1] - skel0[1] + SEPS);
  double drdist = (skel1[2] - skel0[2]) / (dist + SEPS);
  double dy2    = dydx   * dydx;
  double dr2    = drdist * drdist;
  double dydr   = dydx   * drdist;
  double topsq  = sqrt(1. - dr2);
  double denom  = sqrt(1. + dy2);
  //set n0x,n0y,n1x,n1y
  output[0] = (topsq - dydr)          / (denom + SEPS);
  output[1] = (dydx * topsq + drdist) / (denom + SEPS);
  output[2] = (topsq + dydr)          / (denom + SEPS);
  output[3] = (dydx * topsq - drdist) / (denom + SEPS);
}

void reconstruction_s(scalar c, vector n, scalar alpha){
  //Now that we have reconstructed our clean vof cells we need to do reconstruction on our skeleton based cells
  //we can refine our values if needed at the end
  struct skeledata *sd = c.skel;
  struct reconstructTransport *rt = alpha.recont;
  struct reconstructIntermediate *ri = malloc(sizeof(struct reconstructIntermediate));
  scalar fid = c.fid;
  scalar rid = alpha.rid;
  //count allocate & set our intermediate structure
  int ncounts = 0;//count up our neighbors roughly, will give enough space
  int **decoder = malloc(*sd->length * sizeof(int*));
  int *decoderl = calloc(*sd->length , sizeof(int ));
  for(int i = 0; i < *sd->length; i++){
      ncounts += sd->nlengths[i];
      decoderl[i] = sd->nlengths[i];//set our length
      decoder[i]  = calloc(sd->nlengths[i],sizeof(int));//set space for locations 
  }
  double **nri = malloc(ncounts * sizeof(double*));
  int **npair  = malloc(ncounts * sizeof(int*));
  for(int i = 0; i < ncounts; i++){
    nri[i]   = calloc(2*dimension,sizeof(double));
    npair[i] = calloc(2,sizeof(int));
  }
  int *nlength = malloc(sizeof(int));
  *nlength = ncounts;
  ri->normal = nri;
  ri->npairs = npair;
  ri->nlength = nlength;
  ri->decoder = decoder;
  ri->dlength = decoderl;
  //we need to track local index, we can use struct to avoid needing reduction or noauto
#if _MPI
  int *avoidMPI = calloc(1,sizeof(int));
  ri->avoidMPI = avoidMPI;
#else
  int indx = 0;
#endif
  foreach(noauto){
    if(fid[]  > 0.){//Here we compute for each skeleton point
#if _MPI
      int indx = *ri->avoidMPI;
#endif
      int nref = (int)fid[] - 1;//apply our downshift
      double *mainPt = sd->sPts[nref];//grabs our current point
      for(int i = 0; i < sd->nlengths[nref];i++){
        //determine if pair already has been computed, choose is lower val, rchoose is higherval
        int nrefn = sd->neighbors[nref][i];
        int lchoose = nref < nrefn ? nref : nrefn;//holds lower id val
        int uchoose = lchoose == nref ? nrefn : nref;//holds upper idval
        int pass = 1;
        for(int q = 0; q < ri->dlength[nrefn]; q++){
          if(ri->npairs[ri->decoder[nrefn][q]][0] == lchoose &&ri->npairs[ri->decoder[nrefn][q]][1] == uchoose){
            pass = 0;
            break;
          }
        }
        if(pass){
          //loops through each of our neighbors and computes the normal representation of that pair
          double *refPt = sd->sPts[sd->neighbors[nref][i]];//grab our neighboring point
          computeSkelNorm(mainPt,refPt,ri->normal[indx]);
          ri->npairs[indx][0] = lchoose;//keep lower upper ordering
          ri->npairs[indx][1] = uchoose;
          indx++;
        }
      }
#if _MPI
      //if in mpi we need to update our skip/indx value
      *ri->avoidMPI = indx;
#endif
    }
  }
  //now we have the normals for each pair we want to place them inside of our structure n 
  boundary({n,c});
  foreach(noauto){
    if(fid[] > 0.){
      int nref = (int)fid[] - 1;
      int ncount = sd->nlengths[nref];
      double *spt = sd->sPts[nref];
      double *normscompute = malloc(dimension * 2 * sizeof(double));//gather size for both normal averages for this cell
      if(nlengths == 0){
        normNN(point,normscompute,n,c);
      }
      else if(nlengths == 1){
        normON(point,normscompute,n,c,ri,nref);
      }
      else{
        //if more than 1 neighbors we will simply average the 
        normALL(point,normscompute,c,ri,nref);
      }
      //now we have both of our normals averaged from the fields.
#if dimension == 2
      double *bounds = getBounds(x,y,Delta);
#else
      double *bounds = getBounds(x,y,z,Delta);
#endif
      double *n0 = getNormPoint(normscompute,0,spt),*n1 = getNormPoint(normscopmute,1,spt);
      int n0b = pointInsideCell(bounds,n0),n1b = pointInsideCell(bounds,n1b);
#if dimension == 2
      coord n0c = {normscompute[0],normscompute[1]};
      coord n1c = {normscompute[2],normscompute[3]};
#else
      coord n0c = {normscompute[0],normscompute[1],normscompute[2]};
      coord n1c = {normscompute[3],normscompute[4],normscompute[5]};
#endif
      if(n0b && n1b){
        //this means we have 2 interfaces within our cell, and we will set space in our structure
        rid[] = *rt->illength + 1;//upshifts ID for double interface
        rt->interfaceLocations = realloc((*rt->illength + 1 ) * sizeof(int*));//1 size for the over top
        rt->interfacecount = realloc((*rt->illength + 1) * sizeof(int));
        rt->alphas = realloc((*rt->length + 2 ) * sizeof(int));//2 size for under
        rt->normals = realloc((*rt->length + 2) * sizeof(double*));
        rt->normals[*rt->length] = calloc(dimension,sizeof(double));
        rt->normals[*rt->length + 1] = calloc(dimension,sizeof(double));
        rt->interfacecount[*rt->illength] = 2;
        rt->interfaceLocations[*rt->illength] = malloc(2 * sizeof(int));
        rt->interfaceLocations[*rt->illength][0] = *rt->length;
        rt->interfaceLocations[*rt->illength][1] = *rt->length+1;
        //Next add our normal and alpha to structure
        for(int q = 0; q < dimension; q++){
            rt->normals[*rt->length  ][q] = n0[q];
            rt->normals[*rt->length+1][q] = n1[q];
        }
        //conver to coord
        rt->alphas[*rt->length][q] = plane_alpha(c[],n0c);
        rt->alphas[*rt->length][q] = plane_alpha(c[],n1c);
        //finally add to our numbers
        *rt->illength++;
        *rt->length = *rt->length = 2;
      }
      else{
        //if not, we can add a norm if it falls in our cell, if not we add it to the outside placement
        if(n0b){
          //add n0 to our normal
          int q = 0;
          foreach_dimension(){
            n.x[] = normscompute[q++];
          }
          alpha[] = plane_alpha(c[],n0c);
          //and then we set asside n1 to our fill
          ri->ncarry = realloc((*ri->ncarryl+1) * sizeof(double*));
          ri->npoint = realloc((*ri->ncarryl+1) * sizeof(double*));
          ri->ncarry[*ri->ncarryl] = malloc(dimension * sizeof(double));
          ri->npoint[*ri->ncarryl] = malloc(dimension * sizeof(double));
          for(q = 0; q < dimension; q++){
            ri->ncarry[*ri->ncarryl][q] = normcompute[dimension+q];
            ri->npoint[*ri->ncarryl][q] = n1[q];
          }
          *ri->ncarryl++;
        }
        else if(n1b){
          //add n1 to our normal
          int q = 0;
          foreach_dimension(){
            n.x[] = normscompute[dimension+q++];
          }
          alpha[] = plane_alpha(c[],n1c);
          //and then we set asside n1 to our fill
          ri->ncarry = realloc((*ri->ncarryl+1) * sizeof(double*));
          ri->npoint = realloc((*ri->ncarryl+1) * sizeof(double*));
          ri->ncarry[*ri->ncarryl] = malloc(dimension * sizeof(double));
          ri->npoint[*ri->ncarryl] = malloc(dimension * sizeof(double));
          for(q = 0; q < dimension; q++){
            ri->ncarry[*ri->ncarryl][q] = normcompute[q];
            ri->npoint[*ri->ncarryl][q] = n0[q];
          }
          *ri->ncarryl++; 
        }
        else{
          //set both to fill        }
          ri->ncarry = realloc((*ri->ncarryl+2) * sizeof(double*));
          ri->npoint = realloc((*ri->ncarryl+2) * sizeof(double*));
          ri->ncarry[*ri->ncarryl]   = malloc(dimension * sizeof(double));
          ri->ncarry[*ri->ncarryl+1] = malloc(dimension * sizeof(double));
          ri->npoint[*ri->ncarryl]   = malloc(dimension * sizeof(double));
          ri->npoint[*ri->ncarryl+1] = malloc(dimension * sizeof(double));
          for(q = 0; q < dimension; q++){
            ri->ncarry[*ri->ncarryl][q] = normcompute[q];
            ri->ncarry[*ri->ncarryl+1][q] = normcompute[dimension+q];
            ri->npoint[*ri->ncarryl][q] = n0[q];
            ri->npoint[*ri->ncarryl+1][q] = n1[q];
          }
          *ri->ncarryl = *ri->ncarryl + 2; 
      }
    }
  }
  //lastly, addin elements which exist else outside of cell
  for(int i = 0 ; i < *ri-ncarryl;i++){
#if dimension == 2
    Point point = locate(ri->npoint[i][0],ri->npoint[i][1]);
    coord send = {ri->ncarry[i][0],ri->ncarry[i][1]};
#else
    Point point = locate(ri->npoint[i][0],ri->npoint[i][1],ri->npoint[i][2]);
    coord send = {ri->ncarry[i][0],ri->ncarry[i][1],ri->ncarry[i][2]};
#endif
    int q = 0;
    foreach_dimension()
      n.x[] = ri->ncarry[i][q++];
    alpha[] = plane_alpha(c[],send);
  }
  destroyRInter(ri,sd);
#if TREE
  foreach_dimension()
    n.x.refine = n.x.prolongation = refine_injection;
  alpha.n = n;
  alpha.refine = alpha.prolongation = alpha_refine;
#endif
}
