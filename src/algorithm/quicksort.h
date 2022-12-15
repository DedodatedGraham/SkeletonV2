#ifndef QUICK_HEAD
#define QUICK_HEAD
#include "../datastruct/includeds.h"

void swapPoint(struct point *a,struct point *b);
int partition(struct point **points[],int *dim,int left,int right);
void quicksortPoint(struct point **points[],int *dim,int left,int right);


#endif
