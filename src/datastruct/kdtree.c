#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "obj.h"
#include "kdtree.h"


#define MAX 4

struct kdtree{
    struct point *node;//The trees node point
    struct kdtree *left,*right;
    int depth;
    int k;
};

int main(){
    struct kdtree *tree;
    int k=2,d=0,n=1000;
    printf("creating\n");
    tree = kdCreate(&k,&d,&n);
    //Example for loading in points :_)
    FILE *fp;
    fp = fopen("../../../SkeletonV1/SkeleData/Input/interface_points_020000.dat","r");
    int i;
    for(i = 0;i < 2; i++){

    }
    //Then add to tree
    kdLoad(tree);
}
//Creation Methods
struct kdtree *kdCreate(int *k,int *d,int *n){
    //kdCreate will build tree, get splitting dimensions, and propper depth needed for easy loading
    //Will Allocate needed memory, However Will involve alot of Pointers:/
    struct kdtree *tree;
    if(!(tree = malloc(sizeof *tree))){
        return 0;
    }
    tree->depth = *d;
    tree->k = *d % *k;
    if(*n / 2 > MAX){
        int nextsize = ceil((*n - 1) / 2);//take away node point and find max of layer, we maximize size here.
        int nextdep = *d + 1;
        tree->left = kdCreate(k,&nextdep,&nextsize);
        tree->right = kdCreate(k,&nextdep,&nextsize);
    }
    return tree;
}
void kdLoad(struct kdtree *tree,struct point *Pts[]){

}
