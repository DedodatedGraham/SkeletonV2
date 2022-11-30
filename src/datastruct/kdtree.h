#ifndef KDTREE_HEAD
#define KDTREE_HEAD

struct kdtree;

struct kdtree *kdCreate(int *k,int *d,int *n);
void kdLoad(struct kdtree *tree,struct point *Pts[]);

#endif
