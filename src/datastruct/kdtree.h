#ifndef KDTREE_HEAD
#include "obj.h"
#include "../algorithm/quicksort.h"
#define KDTREE_HEAD
struct kdtree;

struct kdtree *kdCreate(struct kdtree *tree,int *k,int *d,int *n);
void kdLoad(struct kdtree tree,struct point **Pts[],int *len);

#endif
