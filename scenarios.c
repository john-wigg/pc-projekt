#include "scenarios.h"

void hotSpot(double *t, int size, int bheight, int bwidth, int imin, int jmin,
             int g) {
  for (int i = g; i < bheight - g; i++) {
    for (int j = g; j < bwidth - g; j++) {
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
void initCustom(double *t, int size, int bheight, int bwidth, int imin,
                int jmin, int g, Conditions cond) {
  for (int i = g; i < bheight - g; i++) {
    for (int j = g; j < bwidth - g; j++) {
      t[i * bwidth + j] = 0.0;
    }
  }
  int mem = (bwidth - 2 * g) * (bheight - 2 * g);
  for (int k = 0; k < mem; k++) {
    int i = cond.idi[k] + g;
    int j = cond.idj[k] + g;
    t[i * bwidth + j] = cond.t[k];
  }
}
