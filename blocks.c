#include "blocks.h"

void updatePoint(double *t1, double *t2, int i, int j, double adjstep,
                 int size) {
  t2[i * size + j] =
      t1[i * size + j] +
      adjstep * (t1[(i + 1) * size + j] + t1[(i - 1) * size + j] +
                 t1[i * size + (j + 1)] + t1[i * size + (j - 1)] -
                 4 * t1[i * size + j]);
}

void swap(double **t1, double **t2) {
  double *temp;
  temp = *t2;
  *t2 = *t1;
  *t1 = temp;
}