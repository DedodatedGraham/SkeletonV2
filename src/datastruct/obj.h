#ifndef OBJECT_HEAD
#define OBJECT_HEAD

//Structure 
struct point {
    double x;
    double y;
    double z;
};

struct skelepoint{
    struct point *pt;
    
    double r;
};

struct intpoint{
    struct point *pt;
    
    double nx;
    double ny;
    double nz;
};

//Creation
extern struct skelepoint* makeSkele2(double x,double y,double r);
extern struct skelepoint* makeSkele3(double x,double y,double z,double r);

extern struct intpoint* makeInt2(double x,double y,double nx,double ny);
extern struct intpoint* makeInt3(double x,double y,double z,double nx,double ny, double nz);

//Conversion
extern struct skelepoint* toPoints(struct intpoint *point,double r);

extern struct intpoint* toPointi2(struct skelepoint *point,double nx, double ny);
extern struct intpoint* toPointi3(struct skelepoint *point,double nx, double ny, double nz);

//destruction
extern void removePoints(struct skelepoint *point); 
extern void removePointi(struct intpoint *point); 

#endif
