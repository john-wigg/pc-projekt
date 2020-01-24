#include "scenarios.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void initTopLeftHotBotRightCold(double *t, int size, int bheight, int bwidth,
                                int imin, int jmin, int g) {
  int i;
  for (i = g; i < bheight - g; i++) {
    int j;
    for (j = g; j < bwidth - g; j++) {
      t[i * bwidth + j] = 0.0;
      if (jmin + j - g == 0) t[i * bwidth + j] = 25.0;
      if (jmin + j - g == size - 1) t[i * bwidth + j] = -25.0;
      if (imin + i - g == 0) t[i * bwidth + j] = 25.0;
      if (imin + i - g == size - 1) t[i * bwidth + j] = -25.0;
    }
  }
}

void initHotSpot(double *t, int size, int bheight, int bwidth, int imin,
                 int jmin, int g) {
  int i;
  for (i = g; i < bheight - g; i++) {
    int j;
    for (j = g; j < bwidth - g; j++) {
      t[i * bwidth + j] = 0.0;
      int distX = (jmin + j - g - size / 2);
      int distY = (imin + i - g - size / 2);
      int rr = distX * distX + distY * distY;
      if (rr < size * size / 100) {
        t[i * bwidth + j] = 25.0;
      }
    }
  }
}

void initRightHot(double *t, int size, int bheight, int bwidth, int imin,
                  int jmin, int g) {
  int i;
  for (i = g; i < bheight - g; i++) {
    int j;
    for (j = g; j < bwidth - g; j++) {
      t[i * bwidth + j] = 0.0;
      if (jmin + j - g == size - 1) t[i * bwidth + j] = 25.0;
    }
  }
}
