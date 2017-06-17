#include <stdio.h>
#include <stdlib.h>
#include "comp.h"

int main(){
  printf("test %f %f\n",fexpf(3), cexpf(3));

  double*x = (double*)malloc(3*sizeof(double));
  x[0] = 0.0;
  x[1] = 1.0;
  x[2] = 2.0;
  printf("%f %f %f\n",x[0], x[1], x[2]);
  x = fexpd_v(x, 3);
  printf("%f %f %f\n",x[0], x[1], x[2]);
  free(x);
  return 0;

}
