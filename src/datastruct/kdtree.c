#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"
//Define for usage of functions


#define MAX 4

struct kdtree{
    struct point *node;//The trees node point
    struct kdtree *left,*right;
    int depth;
    int k;
};

//Creation Methods
struct kdtree *kdCreate(struct kdtree *tree,int *k,int *d,int *n){
    //kdCreate will build tree, get splitting dimensions, and propper depth needed for easy loading
    //Will Allocate needed memory, However Will involve alot of Pointers:/
    if(!(tree = malloc(sizeof *tree))){
        return 0;
    }
    tree->depth = *d;
    tree->k = *d % *k;
    if(*n / 2 > MAX){
        int nextsize = ceil((*n - 1) / 2);//take away node point and find max of layer, we maximize size here.
        int nextdep = *d + 1;
        tree->left = kdCreate(tree->left,k,&nextdep,&nextsize);
        tree->right = kdCreate(tree->right,k,&nextdep,&nextsize);
    }
    return tree;
}

void kdLoad(struct kdtree tree,struct point **Pts[],int *len){
    //Sorts our list
    printf("len passed %d\n",*len);
    int start = 0;
    int stop = *len;
    printf("sorting Points\n");
    printf("Inital p1: x=%lf y=%lf\n",Pts[0]->x,Pts[0].y);
    quicksortPoint(Pts,&tree.k,start,stop);
    printf("After p1: x=%lf y=%lf\n",Pts[0],Pts[0].y);
    
    //Then we choose our nodepoint
    //printf("%lf\n",Pts[len].x);
}
