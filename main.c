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
#include <string.h>

#include "blocks.h"
#include "io.h"
#include "scenarios.h"

int main(int argc, char **argv)
{
  // Initialisiere MPI
  int res = MPI_Init(&argc, &argv);

  if (res != MPI_SUCCESS)
  {
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
  char filename[64];
  char ifilename[64];
  char scenario[64];

  // TODO: beliebige alpha zulassen (wird aber eigentlich nur für das
  // ermitteln der Zeit am Ende benötigt)
  // TODO: beliebige Seitenlängen a zulassen
  // TODO: beliebigen Sicherheitsfaktor für die Schrittweite zulassen
  double alpha;
  double adj;
  double a; // Seitenlänge.

  // Übergabeparameter für [size iter g filename] einlesen
  if (argc != 3 || !sscanf(argv[1], "%s", ifilename) ||
      !sscanf(argv[2], "%s", filename))
  {
    if (rank == 0)
      printf("Nutzung: %s <ifile> <ofile>\n", argv[0]);
    MPI_Finalize();
    return -1;
  }
  readInputFile(&size, &iter, &g, &adj, &alpha, &a, scenario, ifilename);

  double h = a / size;
  double t = 0.0;
  double dt = adj * h * h / 4.0 / alpha;

  // 2 Speicherbereiche für das Wärmefeld
  double *u1, *u2;

  // Anzahl von Blöcken per Zeile und Spalte
  // TODO: Aktuell werden quadratische Blöcke zu Vereinfachung angenommen.
  int num_blocks_per_row;
  int num_blocks_per_col;

  // Finde die optimale Anzahl an Blöcken pro Spalten (sodass die Anzahl Blöcke
  // pro Zeile möglichst groß ist)
  num_blocks_per_row = round(sqrt(num));
  while (num % num_blocks_per_row != 0)
  {
    num_blocks_per_row--;
  }

  num_blocks_per_col = num / num_blocks_per_row;

  if (rank == 0)
  {
    printf("Gitter in %dx%d Blöcke aufgeteilt.\n", num_blocks_per_col,
           num_blocks_per_row);
  }

  // Blockinindizes.
  int bi, bj;
  bi = rank / num_blocks_per_row; // Zeilenindex.
  bj = rank % num_blocks_per_row; // Spaltenindex.

  // Berechne Größe der Blöcke (werden als quadratisch angenommen).
  int bwidth, bheight;
  bheight = ceil((double)size / (double)num_blocks_per_col) + 2 * g;
  bwidth = ceil((double)size / (double)num_blocks_per_row) + 2 * g;

  // Indizes der oberen linken Zelle.
  int imin, jmin;
  imin = (bheight - 2 * g) * bi;
  jmin = (bwidth - 2 * g) * bj;

  int mem = bheight * bwidth * sizeof(double);
  u1 = (double *)malloc(mem);
  u2 = (double *)malloc(mem);

  // Initialisiere das Gitter je nach gewähltem Szenario.
  if (strcmp(scenario, "TOPLEFTHOTBOTRIGHTCOLD") == 0)
  {
    initTopLeftHotBotRightCold(u1, size, bheight, bwidth, imin, jmin, g);
  }
  else if (strcmp(scenario, "HOTSPOT") == 0)
  {
    initHotSpot(u1, size, bheight, bwidth, imin, jmin, g);
  }
  else if (strcmp(scenario, "RIGHTHOT") == 0)
  {
    initRightHot(u1, size, bheight, bwidth, imin, jmin, g);
  }
  else
  {
    initFromImage(u1, size, bheight, bwidth, imin, jmin, g, scenario);
  }

  // Berechne Schwrittweite.
  double adjstep = alpha * adj * 0.25;

  // Iteration im Block.
  MPI_Request req[4];
  MPI_Status stat;
  double *buf;
  buf = (double *)malloc(g * bheight * sizeof(double));
  int l, k, i, j;

  double t1 = MPI_Wtime();
  double time_spent_comm = 0.0;

  for (l = 0; l < iter / g; l++)
  {
    double t1comm = MPI_Wtime();
    // Austauschen der Randbereiche.
    int n_left = bi * num_blocks_per_row + bj - 1;
    int n_right = bi * num_blocks_per_row + bj + 1;
    int n_top = (bi - 1) * num_blocks_per_row + bj;
    int n_bot = (bi + 1) * num_blocks_per_row + bj;

    if (bi > 0)
    {
      // Oberen Rand nach oben senden und unteren Rand von Oben empfangen
      MPI_Isend(u1 + g * bwidth, g * bwidth, MPI_DOUBLE, n_top, T_SENDUP,
                MPI_COMM_WORLD, &req[T_SENDUP]);
      MPI_Recv(u1, g * bwidth, MPI_DOUBLE, n_top, T_SENDDOWN, MPI_COMM_WORLD,
               &stat);
      MPI_Wait(&req[T_SENDUP], &stat);
    }
    if (bi < num_blocks_per_col - 1)
    {
      // Unteren Rand nach unten senden und oberen Rand von unten empfangen
      MPI_Isend(u1 + bwidth * (bheight - 2 * g), g * bwidth, MPI_DOUBLE, n_bot,
                T_SENDDOWN, MPI_COMM_WORLD, &req[T_SENDDOWN]);
      MPI_Recv(u1 + bwidth * (bheight - g), g * bwidth, MPI_DOUBLE, n_bot,
               T_SENDUP, MPI_COMM_WORLD, &stat);
      MPI_Wait(&req[T_SENDDOWN], &stat);
    }

    if (bj > 0)
    {
      // Linken Rand nach links senden und rechten Rand von links empfangen.
      int m, n;
      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          buf[m * g + n] = u1[m * bwidth + g + n];
        }
      }
      MPI_Isend(buf, g * bheight, MPI_DOUBLE, n_left, T_SENDLEFT,
                MPI_COMM_WORLD, &req[T_SENDLEFT]);
      MPI_Recv(buf, g * bheight, MPI_DOUBLE, n_left, T_SENDRIGHT,
               MPI_COMM_WORLD, &stat);
      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          u1[m * bwidth + n] = buf[m * g + n];
        }
      }
      MPI_Wait(&req[T_SENDLEFT], &stat);
    }

    if (bj < num_blocks_per_row - 1)
    {
      // Rechten Rand nach rechts senden und linken Rand von rechts
      int m, n;
      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          buf[m * g + n] = u1[m * bwidth + bwidth - 2 * g + n];
        }
      }
      MPI_Isend(buf, g * bheight, MPI_DOUBLE, n_right, T_SENDRIGHT,
                MPI_COMM_WORLD, &req[T_SENDRIGHT]);
      MPI_Recv(buf, g * bheight, MPI_DOUBLE, n_right, T_SENDLEFT,
               MPI_COMM_WORLD, &stat);
      MPI_Wait(&req[T_SENDRIGHT], &stat);

      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          u1[m * bwidth + bwidth - g + n] = buf[m * g + n];
        }
      }
    }

    double t2comm = MPI_Wtime();
    time_spent_comm += t2comm - t1comm;

    for (k = 0; k < g; k++)
    {
      // Der zu aktualisierende Bereich schrumpft nach jedem Zeitschritt um 1.
      for (i = 1 + k; i < bheight - 1 - k; i++)
      {
        for (j = 1 + k; j < bwidth - 1 - k; j++)
        {
          // Randbereiche mit Dicke g sollen nicht mit aktualisiert werden.
          // TODO: if-Abfragen entfernen
          if (imin + i - g <= 0 || imin + i - g >= size - 1 ||
              jmin + j - g <= 0 || jmin + j - g >= size - 1)
            continue;
          updatePoint(u1, u2, i, j, adjstep, bwidth);
        }
      }
      swap(&u1, &u2);
    }
    t += dt;
  }

  double t2 = MPI_Wtime();
  double time_per_step = (t2 - t1) / iter;
  time_spent_comm /= iter;

  free(buf);

  // Blöcke einsammeln
  double *recvbuf;
  double *u;
  if (rank == 0)
  {
    recvbuf = (double *)malloc(num * bheight * bwidth * sizeof(double));
    u = (double *)malloc(size * size * sizeof(double));
  }

  MPI_Gather(u1, bheight * bwidth, MPI_DOUBLE, recvbuf, bheight * bwidth,
             MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Sortiere das Empfangene Array.
  // TODO: Da findet man sicherlich noch ne schönere Methode
  if (rank == 0)
  {
    int r;
    for (r = 0; r < num; r++)
    {
      int bi, bj;
      bi = r / num_blocks_per_row; // Zeilenindex.
      bj = r % num_blocks_per_row; // Spaltenindex.

      // Indizes der oberen linken Zelle.
      int imin, jmin;
      imin = (bheight - 2 * g) * bi;
      jmin = (bwidth - 2 * g) * bj;

      for (i = g; i < bheight - g; i++)
      {
        if (imin + i - g > size - 1)
          continue; // Falls Blöcke über die Ränder hinausragen.
        for (j = g; j < bwidth - g; j++)
        {
          if (jmin + j - g > size - 1)
            continue;
          u[(imin + i - g) * size + (jmin + j - g)] =
              recvbuf[r * bheight * bwidth + i * bwidth + j];
        }
      }
    }
    printPPMP6(u, size, "out.ppm");
    printf("time per step: %f, time on comm per step: %f, t: %f\n", time_per_step, time_spent_comm, t);
  }

  free(u1);
  free(u2);

  if (rank == 0)
  {
    free(recvbuf);
    free(u);
  }
  MPI_Finalize();
  return 0;
}
