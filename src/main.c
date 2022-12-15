#include "includesrc.c"
int main(){
    printf("Started\n"); 
    struct kdtree tree; 
    
    
    //Example for loading in points :_)
    FILE *fp;
    fp = fopen("../../SkeletonV1/SkeleData/Input/interface_points_020000.dat","r");
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
    //Allocate for loading
    struct intpoint *points[count]; 
    struct point **pts[count]; 
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
        printf("point: %lf %lf %lf %lf\n",x,y,nx,ny);
    }
    //Close
    fclose(fp);
    
    //Then add to tree
    kdCreate(&tree,&k,&d,&count);
    printf("initalized\n");
    i = 0;
    while(i < count){
        printf("adding point %d",i);
        pts[i] = getPointi(points[i]);
        printf(" added: x:%lf y:%lf \n",pts[i].x,pts[i].y);
        i ++;
    }
    printf("converting\n");
    kdLoad(tree,pts,&count);
}
