#include "io.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void printPPM(double *t, int size, const char *filename) {
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

void printPPMP6(double *t, int size, const char *filename) {
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

void printPPMP6MPI(double *t, int size, int bwidth, int bheight, int imin,
                   int jmin, int g, const char *filename, int rank) {
  char header[128];
  sprintf(header, "P6\n%i %i\n255\n", size, size);
  int header_len = strlen(header);
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &fh);

  if (rank == 0) {
    MPI_File_write_at(fh, 0, header, header_len, MPI_BYTE, MPI_STATUS_IGNORE);
  }

  double tmax = 25.0;
  double tmin = -tmax;
  __uint8_t *rgb;
  rgb = (__uint8_t *)malloc((bwidth - 2 * g) * 3 * sizeof(__uint8_t));

  for (int i = g; i < bheight - g; i++) {
    for (int j = g; j < bwidth - g; j++) {
      double val = t[j + i * bwidth];
      int off = 3 * (j - g);  // local offset
      rgb[0 + off] = 0;
      rgb[1 + off] = 0;
      rgb[2 + off] = 0;
      if (val <= tmin) {
        rgb[2 + off] = 255;
      } else if (val >= -25.0 && val < -5) {
        rgb[2 + off] = 255;
        rgb[1 + off] = (__uint8_t)(255 * ((val + 25) / 20));
      } else if (val >= -5 && val <= 0.0) {
        rgb[1 + off] = 255;
        rgb[2 + off] = (__uint8_t)(255 * (1.0 - (val + 5) / 5));
      } else if (val > 0.0 && val <= 5) {
        rgb[1 + off] = 255;
        rgb[0 + off] = (__uint8_t)(255 * ((val) / 5));
      } else if (val > 5 && val < 25.0) {
        rgb[0 + off] = 255;
        rgb[1 + off] = (__uint8_t)(255 * ((25 - val) / 20));
      } else {
        rgb[0 + off] = 255;
      }
    }
    MPI_Offset offset =
        ((imin + i - g) * size + jmin) * 3 * sizeof(__uint8_t) + header_len;
    MPI_File_write_at(fh, offset, rgb, 3 * (bwidth - 2 * g), MPI_UNSIGNED_CHAR,
                      MPI_STATUS_IGNORE);
  }
  MPI_File_close(&fh);
}

void readInputFile(int *size, int *iter, int *g, double *adj, double *alpha,
                   double *a, int *scenario, int *ostep, const char *filename) {
  FILE *f = fopen(filename, "r");

  fscanf(f, "%d\n", size);
  fscanf(f, "%d\n", iter);
  fscanf(f, "%d\n", g);
  fscanf(f, "%lf\n", adj);
  fscanf(f, "%lf\n", alpha);
  fscanf(f, "%lf\n", a);
  fscanf(f, "%d\n", scenario);
  fscanf(f, "%d\n", ostep);

  fclose(f);
}