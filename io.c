#include "io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printPPMP3(double *t, int size, char *filename) {
  FILE *f = fopen(filename, "w");
  fprintf(f, "P3\n%i %i\n255\n", size, size);
  double tmax = 25.0;
  double tmin = -tmax;
  double r, g, b;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      double val = t[j + i * size];
      r = 0;
      g = 0;
      b = 0;
      if (val <= tmin) {
        b = 1.0 * 255.0;
      } else if (val >= -25.0 && val < -5) {
        b = 255 * 1.0;
        g = 255 * ((val + 25) / 20);
      } else if (val >= -5 && val <= 0.0) {
        g = 255 * 1.0;
        b = 255 * (1.0 - (val + 5) / 5);
      } else if (val > 0.0 && val <= 5) {
        g = 255 * 1.0;
        r = 255 * ((val) / 5);
      } else if (val > 5 && val < 25.0) {
        r = 255 * 1.0;
        g = 255 * ((25 - val) / 20);
      } else {
        r = 255 * 1.0;
      }
      fprintf(f, "%i\n%i\n%i\n", (int)r, (int)g, (int)b);
    }
    //      fprintf(f,"\n");
  }
  fclose(f);
}

void printPPMP6(double *t, int size, char *filename) {
  FILE *f = fopen(filename, "w");
  fprintf(f, "P6\n%i %i\n255\n", size, size);
  double tmax = 25.0;
  double tmin = -tmax;
  double r, g, b;
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      double val = t[j + i * size];
      r = 0;
      g = 0;
      b = 0;
      if (val <= tmin) {
        b = 1.0 * 255.0;
      } else if (val >= -25.0 && val < -5) {
        b = 255 * 1.0;
        g = 255 * ((val + 25) / 20);
      } else if (val >= -5 && val <= 0.0) {
        g = 255 * 1.0;
        b = 255 * (1.0 - (val + 5) / 5);
      } else if (val > 0.0 && val <= 5) {
        g = 255 * 1.0;
        r = 255 * ((val) / 5);
      } else if (val > 5 && val < 25.0) {
        r = 255 * 1.0;
        g = 255 * ((25 - val) / 20);
      } else {
        r = 255 * 1.0;
      }
      fprintf(f, "%c%c%c", (char)r, (char)g, (char)b);
    }
    //      fprintf(f,"\n");
  }
  fclose(f);
}

void readGraymap(int **idi, int **idj, double **t, int size, int bheight,
                 int bwidth, int imin, int jmin, int g, char *filename) {
  FILE *f = fopen(filename, "r");
  if (!f) {
    printf("Could not open input file!\n");
    exit(-1);
  }
  char magic[256];
  int iwidth, iheight;
  fscanf(f, "%s\n", magic);
  if (strcmp(magic, "P6") != 0) {
    printf(
        "Nur binärkodieret graustufige PPM-Dateien sind als Input zulässig! "
        "%s\n",
        magic);
    exit(-1);
  }
  fscanf(f, "%d %d\n", &iwidth, &iheight);
  if (iwidth != iheight) {
    printf("Das Bild muss quadratisch sein!");
  }
  int maxval;
  fscanf(f, "%d\n", &maxval);

  // Reserviere die maximal nötige Speichermenge
  // TODO: Weniger Speicher reservieren, wenn weniger gebraucht wird
  int mem = (bwidth - 2 * g) * (bheight - 2 * g);
  *idi = (int *)malloc(mem * sizeof(int));
  *idj = (int *)malloc(mem * sizeof(int));
  *t = (double *)malloc(mem * sizeof(double));

  // TODO: Was, wenn das Bild größer ist als das Gitter?
  int spacing = size / iwidth;
  int i, j;
  int k = 0;
  int val = 0;
  for (i = 0; i < iheight; i++) {
    for (j = 0; j < iwidth; j++) {
      fread(&val, 1, 1, f);
      fread(&val, 1, 1, f);
      fread(&val, 1, 1, f);  // Nur B-Channel wird gelesen
      if (i * spacing >= imin && (i + 1) * spacing < imin + bheight - 2 * g) {
        if (j * spacing >= jmin && (j + 1) * spacing < jmin + bwidth - 2 * g) {
          if (val != 0) {
            for (int il = 0; il < spacing; il++) {
              for (int jl = 0; jl < spacing; jl++) {
                (*idi)[k] = i * spacing + il - imin;
                (*idj)[k] = j * spacing + jl - jmin;
                (*t)[k] = (double)val;
                k++;
              }
            }
          }
        }
      }
    }
  }
}