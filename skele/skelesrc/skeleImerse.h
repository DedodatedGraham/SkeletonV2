#ifndef _skeleIm_
coord skeleton = {0.,0.};
extern double Re;
extern double Cd;
extern double D;
extern double u0;

double computeDis(coord skele, coord location,double r){
  double ret = 0.;
  foreach_dimension()
      ret+=sq(skele.x - location.x);
  ret = sqrt(ret);
  return ret <= r ? ret : 0.; 
}
double computeGausian(double sigma, double xloc){
  return exp(-1 * sq(xloc) / (2*sq(sigma))) * (1 / (sigma * sqrt(2.*3.1415)));
}
void skeleImerseForce(Point point,face vector av, coord force, coord skele,double r){
  if(is_refined(cell)){
    //cell has children so we go deeper
    foreach_child()
      skeleImerseForce(point,av,force,skele,r);
  }
  else{
    //When we have reached our point we want to apply our force
    double sig = r/3;
    coord loc = {x , y};
    double nx =(x-skele.x),ny=(y-skele.y);
    double sqs=sqrt(sq(nx)+sq(ny));
    nx=nx/sqs,ny=ny/sqs;
    double cdis = computeDis(skele,loc,r);//(x - skele.x);
    if(cdis){
        //add in xforces
        double gauval = computeGausian(sig,cdis); 
        av.x[0] += force.x*gauval*nx;//apply gausian drag x-
        av.x[1] += force.x*gauval*nx;//apply gausian drag x+
        av.y[0] += force.x*gauval*ny;//apply gausian drag y-
        av.y[1] += force.x*gauval*ny;//apply gausian drag y+
    }
  }
}
void scopeImerse(Point point,int lnow,coord skele, double r, face vector av){
  int pass = 1;
  foreach_child(){
    int lpass = 1;
    foreach_dimension(){
      double lb = x - Delta/2.;
      double rb = x + Delta/2.;
      if(skele.x < lb || skele.x > rb)lpass = -1;
      else if(lpass != -1 && skele.x - r > lb && skele.x + r < rb)lpass = 0;
    }
    if(!lpass){
      scopeImerse(point,lnow+1,skele,r,av);//only go when 0 :)
      pass = 0;
    }
  }
  //when pass we apply forces into cells below 
  if(pass>0){
    coord force = {Cd*0.5*(3.1415*sq(r)), 0.};
    skeleImerseForce(point,av,force,skele,r);
  }
}
event acceleration(i++){
  printf("accel @ i = %d\n",i);
  //when acceleating skeleton point we first want to locate relevant regions of our skeleton
  double r = D/2.;//ensure r
  //we run a short search for the wanted level starting from our top level 
  face vector av = a;
  foreach_level(0){
    scopeImerse(point,0,skeleton,r,av);
  }
}
#endif
