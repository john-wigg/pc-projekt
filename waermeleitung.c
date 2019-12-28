#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Initialisiere Wärmefeld mit Startwerten:
// innen: 0.0
// Rand:
// links/oben warm=25.0
// rechts/unten kalt=-25.0
void init(double *t, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      t[j + i * size] = 0.0;
      if (j == 0) t[i * size] = 25.0;
      if (j == size - 1) t[j + i * size] = -25.0;
      if (i == 0) t[j + i * size] = 25.0;
      if (i == size - 1) t[j + i * size] = -25.0;
    }
  }
}

// Euler-Vorwärts für einen einzelnen Gitterpunkt.
// adjstep ist die Schrittweite des Euler-Verfahrens.
void updateCell(double *t1, double *t2, int i, int j, double adjstep,
                int size) {
  t2[i * size + j] =
      t1[i * size + j] +
      adjstep * (t1[(i + 1) * size + j] + t1[(i - 1) * size + j] +
                 t1[i * size + (j + 1)] + t1[i * size + (j - 1)] -
                 4 * t1[i * size + j]);
}

// Vertausche zwei double-Pointer.
// Wird benutzt, um das aktualisierte und das alte Gitter auszutauschen.
void swap(double **t1, double **t2) {
  double *temp;
  temp = *t2;
  *t2 = *t1;
  *t1 = temp;
}

// Ausgabe des Feldes t als PPM (Portable Pix Map) in filename
// mit schönen Farben
void printResult(double *t, int size, char *filename) {
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
        g = 255 * (1.0 - (val - 25) / 20);
      } else {
        r = 255 * 1.0;
      }
      fprintf(f, "%i\n%i\n%i\n", (int)r, (int)g, (int)b);
    }
    //      fprintf(f,"\n");
  }
  fclose(f);
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  // Größe des Feldes
  int size;
  // Anzahl Iterationen
  int iter;
  // Geisterzonenbreite
  int g;
  // Ausgabedatei
  char filename[256];

  // Übergabeparameter für [size iter g filename] einlesen
  if (argc != 5) {
    printf("Nutzung: %s <size> <iter> <g> <filename>\n", argv[0]);
    MPI_Finalize();
    return -1;
  }

  if (!sscanf(argv[1], "%d", &size) || !sscanf(argv[2], "%d", &iter) ||
      !sscanf(argv[3], "%d", &g) || !sscanf(argv[4], "%s", filename)) {
    printf("Nutzung: %s <size> <iter> <g> <filename>\n", argv[0]);
    MPI_Finalize();
    return -1;
  }

  // 2 Speicherbereiche für das Wärmefeld
  double *u1, *u2;

  // Größe des Speicherbereiches
  int mem = size * size * sizeof(double);
  // Allokiere Speicher auf Host
  u1 = (double *)malloc(mem);
  u2 = (double *)malloc(mem);
  // Initialisiere Speicher
  init(u1, size);

  // TODO: Implementieren Sie ein paralleles Programm,
  //      das die Temperaturen u_k mit einer beliebigen
  //      Anzahl von p Prozessoren unter Verwendung von MPI berechnet.
  //      Der Austausch der Randbereiche soll alle g Schritte passieren.

  // Gib das Ergebnis aus
  printResult(u1, size, filename);

  MPI_Finalize();
  return 0;
}
