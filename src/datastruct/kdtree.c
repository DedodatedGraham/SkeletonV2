#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "obj.c"
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
    //Example for loading in points :_)
    FILE *fp;
    fp = fopen("../../../SkeletonV1/SkeleData/Input/interface_points_020000.dat","r");
    //Check    
    if (fp == NULL){
        return 0;
    }
    //CountLines
    int count = 0;
    char c;
    for(c = getc(fp);c != EOF; c = getc(fp)){
        if(c == '\n'){
            count++;
        }
    }
    printf("count: %d\n",count);
    //Allocate for loading
    struct intpoint *points[count]; 
    //Loadin
    fseek(fp, 0, SEEK_SET);
    int k = 2;
    int d = 0;
    printf("loading\n");
    int i;
    for(i = 0;i<count;i++){
        //For 2D
        double x,y,nx,ny;
        fscanf(fp,"%lf %lf %lf %lf",&x,&y,&nx,&ny);
        points[i] = makeInt2(x,y,nx,ny); 
    }
    //Close
    fclose(fp);
    //Then add to tree
    printf("creating\n");
    tree = kdCreate(&k,&d,&count);
    printf("Len points: %lu\n",sizeof(points)/sizeof(struct intpoint));
    //kdLoad(tree);
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
