#ifndef OBJECT_HEAD
#define OBJECT_HEAD

//Structure
struct point;

struct skelepoint;

struct intpoint;

//Creation
struct skelepoint* makeSkele2(double x,double y,double r);
struct skelepoint* makeSkele3(double x,double y,double z,double r);

struct intpoint* makeInt2(double x,double y,double nx,double ny);
struct intpoint* makeInt3(double x,double y,double z,double nx,double ny, double nz);

//Conversion
struct skelepoint* toPoints(struct intpoint *point,double r);

struct intpoint* toPointi2(struct skelepoint *point,double nx, double ny);
struct intpoint* toPointi3(struct skelepoint *point,double nx, double ny, double nz);

//destruction
void removePoints(struct skelepoint *point); 
void removePointi(struct intpoint *point); 

#endif
