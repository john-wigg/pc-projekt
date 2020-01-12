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

/**
 * @brief Enthält die Konfiguration für die Simulation
 */
typedef struct config_s
{
  double alpha; //!< thermische Diffusitivität
  double adj;   //!< Sicherheitsfaktor für den Zeitschritt (<1.0 für stabile Simulation)
  double a;     //<! Seitenlänge des Gitters
  int size;     //<! Zellen pro Gitterseite
  int iter;     //<! Maximale Iterationsschritte
  int g;        //<! Breite der Geisterzonen in Zellen
  int scenario; //<! Gewähltes Szeneria, d.h. Anfangs- und Randbedingungen
} t_config;

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

  // Ausgabedatei
  char filename[64];
  char ifilename[64];

  // Es könnte zwar jeder Prozess die Inputdatei einlesen, um Leseoperationen zu sparen,
  // lassen wir jedoch nur Rang 0 die Konfiguration einlesen und verteilen diese dann
  // per MPI_Bcast auf alle anderen Prozesse.

  t_config config; // Jeder Prozess definiert ein Struct, welches die Konfiguration enthält.

  // Serialisieren des Structs mittel MPI_Type_create_struct.
  const int count = 7;                                  // Anzahl der Blöcke im Struct hier entspricht jeder Block einem Element.
  int array_of_blocklengths[7] = {1, 1, 1, 1, 1, 1, 1}; // Längen der Blöcke im Struct.
  MPI_Datatype array_of_types[7] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT,
                                    MPI_INT, MPI_INT, MPI_INT}; // Typen der Blöcke im Struct.
  MPI_Aint array_of_displacements[count];                       // Offsets der Blöcke im Struct.
  array_of_displacements[0] = offsetof(t_config, alpha);
  array_of_displacements[1] = offsetof(t_config, adj);
  array_of_displacements[2] = offsetof(t_config, a);
  array_of_displacements[3] = offsetof(t_config, size);
  array_of_displacements[4] = offsetof(t_config, iter);
  array_of_displacements[5] = offsetof(t_config, g);
  array_of_displacements[6] = offsetof(t_config, scenario);
  MPI_Datatype config_type; // Definition des serialisierten Typen
  MPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements,
                         array_of_types, &config_type);
  MPI_Type_commit(&config_type); // Typen committen.

  // Übergabeparameter für [size iter g filename] einlesen
  if (argc != 3 || !sscanf(argv[1], "%s", ifilename) ||
      !sscanf(argv[2], "%s", filename))
  {
    if (rank == 0)
      printf("Nutzung: %s <ifile> <ofile>\n", argv[0]);
    MPI_Finalize();
    return -1;
  }

  // Rand 0 liest die Inputdatei und veteilt die Konfiguration an alle Prozesse.
  if (rank == 0)
  {
    readInputFile(&config.size, &config.iter, &config.g, &config.adj,
                  &config.alpha, &config.a, &config.scenario, ifilename);
  }
  MPI_Bcast(&config, 1, config_type, 0, MPI_COMM_WORLD);

  MPI_Type_free(&config_type); // Typen löschen.

  int size;
  int iter;
  int g;
  double alpha;
  double adj;
  double a;
  int scenario;

  size = config.size;
  iter = config.iter;
  g = config.g;
  alpha = config.alpha;
  adj = config.adj;
  a = config.a;
  scenario = config.scenario;

  double h = a / size;                   // Zellenbreite
  double dt = adj * h * h / 4.0 / alpha; // Zeitschritt
  double t = 0.0;                        // Zeit

  // Speicherbereiche für das 2D-Gitter des Blockes.
  double *u1, *u2;

  // Anzahl von Blöcken (Prozessen) per Zeile und Spalte.
  int num_blocks_per_row;
  int num_blocks_per_col;

  // Finde die optimale Anzahl an Blöcken pro Zeile/Spalte:
  // Die Differenz zwischen num_blocks_per_row und num_blocks_per_col sollte möglichst klein sein.
  // num_blocks_per_col sollte größer als num_blocks_per_row sein, da das die Übertragung der Geisterzonen
  // einfacher macht.
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

  // Berechne Größe des Blocks (inklusive Geisterzellen).
  int bwidth, bheight;
  bheight = ceil((double)size / (double)num_blocks_per_col) + 2 * g;
  bwidth = ceil((double)size / (double)num_blocks_per_row) + 2 * g;

  // Indizes der oberen linken Zelle im Block.
  int imin, jmin;
  imin = (bheight - 2 * g) * bi;
  jmin = (bwidth - 2 * g) * bj;

  // Reserviere des Speicher der Zelle.
  int mem = bheight * bwidth * sizeof(double);
  u1 = (double *)malloc(mem);
  u2 = (double *)malloc(mem);

  // Initialisiere das Gitter je nach gewähltem Szenario.
  if (scenario == 0)
  {
    initTopLeftHotBotRightCold(u1, size, bheight, bwidth, imin, jmin, g);
  }
  else if (scenario == 1)
  {
    initHotSpot(u1, size, bheight, bwidth, imin, jmin, g);
  }
  else if (scenario == 2)
  {
    initRightHot(u1, size, bheight, bwidth, imin, jmin, g);
  }
  /*
  else {
    initFromImage(u1, size, bheight, bwidth, imin, jmin, g, scenario);
  }
  */

  // Berechne Schwrittweite.
  // adj ist dabei ein Sicherheitsfaktor, der für ein stabiles Verfahren <= 1.0 sein muss.
  double adjstep = alpha * adj * 0.25;

  // Iteration im Block.
  MPI_Request req[4];
  MPI_Request req_recv[4];
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
    // TODO: Ecken austauschen
    int n_left = bi * num_blocks_per_row + bj - 1;
    int n_right = bi * num_blocks_per_row + bj + 1;
    int n_top = (bi - 1) * num_blocks_per_row + bj;
    int n_bot = (bi + 1) * num_blocks_per_row + bj;

    if (bi > 0)
    {
      // Oberen Rand nach oben senden und unteren Rand von Oben empfangen
      MPI_Issend(u1 + g * bwidth, g * bwidth, MPI_DOUBLE, n_top, T_SENDUP,
                 MPI_COMM_WORLD, &req[T_SENDUP]);
      MPI_Irecv(u1, g * bwidth, MPI_DOUBLE, n_top, T_SENDDOWN, MPI_COMM_WORLD,
                &req_recv[T_SENDUP]);
    }
    if (bi < num_blocks_per_col - 1)
    {
      // Unteren Rand nach unten senden und oberen Rand von unten empfangen
      MPI_Issend(u1 + bwidth * (bheight - 2 * g), g * bwidth, MPI_DOUBLE, n_bot,
                 T_SENDDOWN, MPI_COMM_WORLD, &req[T_SENDDOWN]);
      MPI_Irecv(u1 + bwidth * (bheight - g), g * bwidth, MPI_DOUBLE, n_bot,
                T_SENDUP, MPI_COMM_WORLD, &req_recv[T_SENDDOWN]);
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
      MPI_Issend(buf, g * bheight, MPI_DOUBLE, n_left, T_SENDLEFT,
                 MPI_COMM_WORLD, &req[T_SENDLEFT]);
      MPI_Irecv(buf, g * bheight, MPI_DOUBLE, n_left, T_SENDRIGHT,
                MPI_COMM_WORLD, &req_recv[T_SENDLEFT]);
      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          u1[m * bwidth + n] = buf[m * g + n];
        }
      }
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
      MPI_Issend(buf, g * bheight, MPI_DOUBLE, n_right, T_SENDRIGHT,
                 MPI_COMM_WORLD, &req[T_SENDRIGHT]);
      MPI_Irecv(buf, g * bheight, MPI_DOUBLE, n_right, T_SENDLEFT,
                MPI_COMM_WORLD, &req_recv[T_SENDRIGHT]);

      for (m = 0; m < bheight; m++)
      {
        for (n = 0; n < g; n++)
        {
          u1[m * bwidth + bwidth - g + n] = buf[m * g + n];
        }
      }
    }
    if (bi > 0)
    {
      MPI_Wait(&req[T_SENDUP], &stat);
      MPI_Wait(&req_recv[T_SENDUP], &stat);
    }
    if (bi < num_blocks_per_col - 1)
    {
      MPI_Wait(&req[T_SENDDOWN], &stat);
      MPI_Wait(&req_recv[T_SENDDOWN], &stat);
    }
    if (bj > 0)
    {
      MPI_Wait(&req[T_SENDLEFT], &stat);
      MPI_Wait(&req_recv[T_SENDLEFT], &stat);
    }
    if (bj < num_blocks_per_row - 1)
    {
      MPI_Wait(&req[T_SENDRIGHT], &stat);
      MPI_Wait(&req_recv[T_SENDRIGHT], &stat);
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
    printf("time per step: %f, time on comm per step: %f, t: %f\n",
           time_per_step, time_spent_comm, t);
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
