#include "blocks.h"

void initBlock(double *t, int size, int bheight, int bwidth, int imin, int jmin,
               int g) {
  for (int i = g; i < bheight - g; i++) {
    for (int j = g; j < bwidth - g; j++) {
      t[i * bwidth + j] = 0.0;
      if (jmin + j - g == 0) t[i * bwidth + j] = 25.0;
      if (jmin + j - g == size - 1) t[i * bwidth + j] = -25.0;
      if (imin + i - g == 0) t[i * bwidth + j] = 25.0;
      if (imin + i - g == size - 1) t[i * bwidth + j] = -25.0;
    }
  }
}

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