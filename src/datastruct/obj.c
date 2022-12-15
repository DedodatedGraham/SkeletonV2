#include <stdio.h>
#include <stdlib.h>
#include "obj.h"


//Point Creation
struct skelepoint** makeSkele2(double x,double y,double r){
    //define points
    struct point *pt;
    struct skelepoint *spt;
    //allocate memory
    if(!(pt = malloc(sizeof pt))){
        return 0;
    }
    if(!(spt = malloc(sizeof spt))){
        return 0;
    }
    //assign variables
    pt->x = x;
    pt->y = y;
    spt->pt = pt;
    spt->r = r;
    return spt;
}
struct skelepoint** makeSkele3(double x,double y,double z,double r){
    //define points
    struct point *pt;
    struct skelepoint *spt;
    //allocate memory
    if(!(pt = malloc(sizeof pt))){
        return 0;
    }
    if(!(spt = malloc(sizeof spt))){
        return 0;
    }
    //assign variables
    pt->x = x;
    pt->y = y;
    pt->z = z;
    spt->pt = pt;
    spt->r = r;
    return spt;
} 
struct intpoint** makeInt2(double x,double y,double nx,double ny){
    //define points
    struct point *pt;
    struct intpoint *ipt;
    //allocate memory
    if(!(pt = malloc(sizeof pt))){
        return 0;
    }
    if(!(ipt = malloc(sizeof ipt))){
        return 0;
    }
    //assign variables
    pt->x = x;
    pt->y = y;
    ipt->pt = pt;
    ipt->nx = nx;
    ipt->ny = ny;
    return ipt;
}
struct intpoint** makeInt3(double x,double y,double z,double nx, double ny,double nz){
    //define points
    struct point *pt;
    struct intpoint *ipt;
    //allocate memory
    if(!(pt = malloc(sizeof pt))){
        return 0;
    }
    if(!(ipt = malloc(sizeof ipt))){
        return 0;
    }
    //assign variables
    pt->x = x;
    pt->y = y;
    pt->z = z;
    ipt->pt = pt;
    ipt->nx = nx;
    ipt->ny = ny;
    ipt->nz = nz;
    return ipt;
}

//Simple Conversion
//from one type of point to another, for when done with the normals and wanting to save a radius
struct skelepoint* toPoints(struct intpoint *point,double r){
    struct skelepoint *spt;
    if(!(spt = malloc(sizeof spt))){
        return 0;
    }
    spt->pt = point->pt;
    spt->r = r;
    return spt;
}
struct intpoint* toPointi2(struct skelepoint *point,double nx,double ny){
    struct intpoint *ipt;
    if(!(ipt = malloc(sizeof ipt))){
        return 0;
    }
    ipt->pt = point->pt;
    ipt->nx = nx;
    ipt->ny = ny;
    return ipt;
}

struct point getPointi(struct intpoint *point){
    struct point holdpoint; 
    holdpoint = *point->pt;
    return holdpoint;
}


//For Point Deletion
void removePoints(struct skelepoint *point){
    if(point){
        free(point);        
    }
}
void removePointi(struct intpoint *point){
    if(point){
        free(point);        
    }
}

