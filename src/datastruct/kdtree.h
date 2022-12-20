#ifndef KDTREE_HEAD
#define KDTREE_HEAD
#include <math.h>
#include "../algorithm/includea.h"
#include "obj.h"
struct kdtree;
//Creation Methods
struct kdtree *kdCreate(struct kdtree *tree,int *k,int *d,int *n);
void kdLoad(struct kdtree tree,struct point **Pts[],int *len);
#endif
