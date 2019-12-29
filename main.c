/**
 * @file main.c
 * @author Wigg
 * @date 29 December 2019
 * @brief Main-Datei des Projektes.
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "blocks.h"
#include "io.h"

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
  // * Die globalen Ränder dürfen nicht mit aktulaisiert werden.
  MPI_Request req;
  MPI_Status stat;
  double *buf;
  buf = (double *)malloc(g * bsize * sizeof(double));
  int l, k, i, j;
  for (l = 0; l < iter / g; l++) {
    // TODO: Randbereiche austauschen.
    int n_left = bi * num_blocks_per_row + bj - 1;
    int n_right = bi * num_blocks_per_row + bj + 1;
    int n_top = (bi - 1) * num_blocks_per_row + bj;
    int n_bot = (bi + 1) * num_blocks_per_row + bj;
    if (bi > 0) {
      // Oberen Rand nach oben senden und unteren Rand von Oben empfangen
      MPI_Isend(u1 + g * bsize, g * bsize, MPI_DOUBLE, n_top, T_SENDUP,
                MPI_COMM_WORLD, &req);
      MPI_Recv(u1, g * bsize, MPI_DOUBLE, n_top, T_SENDDOWN, MPI_COMM_WORLD,
               &stat);
    }
    if (bi < num_blocks_per_col - 1) {
      // Unteren Rand nach unten senden und oberen Rand von unten empfangen
      MPI_Isend(u1 + bsize * (bsize - 2 * g), g * bsize, MPI_DOUBLE, n_bot,
                T_SENDDOWN, MPI_COMM_WORLD, &req);
      MPI_Recv(u1 + bsize * (bsize - g), g * bsize, MPI_DOUBLE, n_bot, T_SENDUP,
               MPI_COMM_WORLD, &stat);
    }

    if (bj > 0) {
      // Linken Rand nach links senden und rechten Rand von links empfangen.
      int m, n;
      for (m = 0; m < bsize; m++) {
        for (n = 0; n < g; n++) {
          buf[m * g + n] = u1[m * bsize + g + n];
        }
      }
      MPI_Isend(buf, g * bsize, MPI_DOUBLE, n_left, T_SENDLEFT, MPI_COMM_WORLD,
                &req);
      MPI_Recv(buf, g * bsize, MPI_DOUBLE, n_left, T_SENDRIGHT, MPI_COMM_WORLD,
               &stat);
      for (m = 0; m < bsize; m++) {
        for (n = 0; n < g; n++) {
          u1[m * bsize + n] = buf[m * g + n];
        }
      }
    }

    if (bj < num_blocks_per_row - 1) {
      // Rechten Rand nach rechts senden und linken Rand von rechts
      int m, n;
      for (m = 0; m < bsize; m++) {
        for (n = 0; n < g; n++) {
          buf[m * g + n] = u1[m * bsize + bsize - 2 * g + n];
        }
      }
      MPI_Isend(buf, g * bsize, MPI_DOUBLE, n_right, T_SENDRIGHT,
                MPI_COMM_WORLD, &req);
      MPI_Recv(buf, g * bsize, MPI_DOUBLE, n_right, T_SENDLEFT, MPI_COMM_WORLD,
               &stat);

      for (m = 0; m < bsize; m++) {
        for (n = 0; n < g; n++) {
          u1[m * bsize + bsize - g + n] = buf[m * g + n];
        }
      }
    }

    for (k = 0; k < g; k++) {
      for (i = 1; i < bsize - 1 - k; i++) {
        for (j = 1; j < bsize - 1 - k; j++) {
          // Randbereiche mit Dicke g sollen nicht mit aktualisiert werden.
          if (imin + i - g == 0 || imin + i - g >= size - 1 ||
              jmin + j - g == 0 || jmin + j - g >= size - 1)
            continue;
          updatePoint(u1, u2, i, j, adjstep, bsize);
        }
      }
      swap(&u1, &u2);
    }
  }

  free(buf);

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

  free(u1);
  free(u2);

  if (rank == 0) {
    free(recvbuf);
    free(u);
  }
  MPI_Finalize();
  return 0;
}
