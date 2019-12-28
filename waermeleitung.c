#include <math.h>
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

// Initialisiere einzelne Blöcke des Wärmefeldes mit Startwerten.
// innen: 0.0
// Rand:
// links/oben warm=25.0
// rechts/unten kalt=-25.0
// Größe der Blöcke bsize sowie Indizes der oberen linken Zelle imin, jmin
// werden benötigt.
void initBlock(double *t, int size, int bsize, int imin, int jmin, int g) {
  for (int i = g; i < bsize - g; i++) {
    for (int j = g; j < bsize - g; j++) {
      t[i * bsize + j] = 0.0;
      if (jmin + j - g == 0) t[i * bsize + j] = 25.0;
      if (jmin + j - g == size - 1) t[i * bsize + j] = -25.0;
      if (imin + i - g == 0) t[i * bsize + j] = 25.0;
      if (imin + i - g == size - 1) t[i * bsize + j] = -25.0;
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

int main(int argc, char **argv) {
  // Initialisiere MPI
  int res = MPI_Init(&argc, &argv);

  if (res != MPI_SUCCESS) {
    MPI_Abort(MPI_COMM_WORLD, res);
  }

  // Größe und Rang in World
  int num, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &num);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Größe des Feldes
  int size;
  // Anzahl Iterationen
  int iter;
  // Geisterzonenbreite
  int g;
  // Ausgabedatei
  char filename[256];

  // TODO: beliebige alpha zulassen (wird aber eigentlich nur für das
  // ermitteln der Zeit am Ende benötigt)
  // TODO: beliebige Seitenlängen a zulassen
  // TODO: beliebigen Sicherheitsfaktor für die Schrittweite zulassen
  double alpha = 1.0;
  double a = 1.0;
  double adj = 1.0;

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

  // TODO: Implementieren Sie ein paralleles Programm,
  //      das die Temperaturen u_k mit einer beliebigen
  //      Anzahl von p Prozessoren unter Verwendung von MPI berechnet.
  //      Der Austausch der Randbereiche soll alle g Schritte passieren.

  // Anzahl von Blöcken per Zeile und Spalte
  // TODO: Aktuell werden quadratische Blöcke zu Vereinfachung angenommen.
  int num_blocks_per_row;
  int num_blocks_per_col;

  num_blocks_per_row = (int)sqrt(num);
  num_blocks_per_col = num_blocks_per_row;

  // Blockinindizes.
  int bi, bj;
  bi = rank / num_blocks_per_col;  // Zeilenindex.
  bj = rank % num_blocks_per_row;  // Spaltenindex.

  // Berechne Größe der Blöcke (werden als quadratisch angenommen).
  int bsize;
  bsize = ceil((double)size / (double)num_blocks_per_col) + 2 * g;

  // Indizes der oberen linken Zelle.
  int imin, jmin;
  imin = (bsize - 2 * g) * bi;
  jmin = (bsize - 2 * g) * bj;

  printf("-------------------------------\n");
  printf("Block %d :: %d | %d :: %d | %d :: %d\n", rank, bi, bj, imin, jmin,
         bsize);

  int mem = bsize * bsize * sizeof(double);
  u1 = (double *)malloc(mem);
  u2 = (double *)malloc(mem);

  initBlock(u1, size, bsize, imin, jmin, g);

  // Berechne Schwrittweite.
  double adjstep = adj * 0.25;

  // Iteration im Block.
  // * Nach jedem Schritt wird der zu aktualisierende Block um 1 verkleinert.
  // * TODO: Nach g Schritten müssen die Randbereiche ausgetauscht werden.
  // * TODO: Die globalen Ränder dürfen nicht mit aktulaisiert werden.
  int l, k, i, j;
  for (l = 0; l < iter / g; l++) {
    // TODO: Randbereiche austauschen.
    for (k = 0; k < g; k++) {
      for (i = 1; i < bsize - 1 - k; i++) {
        for (j = 1; j < bsize - 1 - k; j++) {
          // Randbereiche mit Dicke g sollen nicht mit aktualisiert werden.
          if (imin + i - g == 0 || imin + i - g == size - 1 ||
              jmin + j - g == 0 || jmin + j - g == size - 1)
            continue;
          updateCell(u1, u2, i, j, adjstep, bsize - k);
        }
      }
      swap(&u1, &u2);
    }
  }

  // Blöcke einsammeln
  double *recvbuf;
  double *u;
  if (rank == 0) {
    recvbuf = (double *)malloc(num * bsize * bsize * sizeof(double));
    u = (double *)malloc(size * size * sizeof(double));
  }

  MPI_Gather(u1, bsize * bsize, MPI_DOUBLE, recvbuf, bsize * bsize, MPI_DOUBLE,
             0, MPI_COMM_WORLD);

  // Sortiere das Empfangene Array.
  // TODO: Da findet man sicherlich noch ne schönere Methode
  if (rank == 0) {
    int r;
    for (r = 0; r < num; r++) {
      int bi, bj;
      bi = r / num_blocks_per_col;  // Zeilenindex.
      bj = r % num_blocks_per_row;  // Spaltenindex.

      // Indizes der oberen linken Zelle.
      int imin, jmin;
      imin = (bsize - 2 * g) * bi;
      jmin = (bsize - 2 * g) * bj;

      for (i = g; i < bsize - g; i++) {
        if (imin + i - g > size - 1)
          continue;  // Falls Blöcke über die Ränder hinausragen.
        for (j = g; j < bsize - g; j++) {
          if (jmin + j - g > size - 1) continue;
          u[(imin + i - g) * size + (jmin + j - g)] =
              recvbuf[r * bsize * bsize + i * bsize + j];
        }
      }
    }
    printResult(u, size, "out.ppm");
  }

  // Gib das Ergebnis aus
  char bfname[256];
  sprintf(bfname, "out%d.ppm", rank);
  printResult(u1, bsize, bfname);

  MPI_Finalize();
  return 0;
}
