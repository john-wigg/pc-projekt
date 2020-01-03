#ifndef SCENARIOS_H_
#define SCENARIOS_H_

void hotSpot(double *t, int size, int bheight, int bwidth, int imin, int jmin,
             int g);

typedef struct {
  int *idi;
  int *idj;
  double *t;
} Conditions;

void initCustom(double *t, int size, int bheight, int bwidth, int imin,
                int jmin, int g, Conditions cond);

#endif  // SCENARIOS_H_