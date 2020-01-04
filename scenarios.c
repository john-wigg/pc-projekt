#include "scenarios.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void initTopLeftHotBotRightCold(double *t, int size, int bheight, int bwidth,
                                int imin, int jmin, int g)
{
  for (int i = g; i < bheight - g; i++)
  {
    for (int j = g; j < bwidth - g; j++)
    {
      t[i * bwidth + j] = 0.0;
      if (jmin + j - g == 0)
        t[i * bwidth + j] = 25.0;
      if (jmin + j - g == size - 1)
        t[i * bwidth + j] = -25.0;
      if (imin + i - g == 0)
        t[i * bwidth + j] = 25.0;
      if (imin + i - g == size - 1)
        t[i * bwidth + j] = -25.0;
    }
  }
}

void initHotSpot(double *t, int size, int bheight, int bwidth, int imin,
                 int jmin, int g)
{
  for (int i = g; i < bheight - g; i++)
  {
    for (int j = g; j < bwidth - g; j++)
    {
      t[i * bwidth + j] = 0.0;
      int distX = (jmin + j - g - size / 2);
      int distY = (imin + i - g - size / 2);
      int rr = distX * distX + distY * distY;
      if (rr < size * size / 100)
      {
        t[i * bwidth + j] = 25.0;
      }
    }
  }
}

void initRightHot(double *t, int size, int bheight, int bwidth, int imin,
                  int jmin, int g)
{
  for (int i = g; i < bheight - g; i++)
  {
    for (int j = g; j < bwidth - g; j++)
    {
      t[i * bwidth + j] = 0.0;
      if (jmin + j - g == size - 1)
        t[i * bwidth + j] = 25.0;
    }
  }
}

void initFromImage(double *t, int size, int bheight, int bwidth, int imin,
                   int jmin, int g, char *filename)
{
  FILE *f = fopen(filename, "r");
  if (!f)
  {
    printf("Datei konnte nicht geöffnet werden!\n");
    exit(-1);
  }
  char magic[256];
  int iwidth, iheight;
  fscanf(f, "%s\n", magic);
  if (strcmp(magic, "P6") != 0)
  {
    printf(
        "Nur binärkodierte PPM-Dateien sind als Input zulässig!\n");
    exit(-1);
  }
  fscanf(f, "%d %d\n", &iwidth, &iheight);
  if (iwidth != iheight)
  {
    printf("Das Bild muss quadratisch sein!\n");
  }
  int maxval;
  fscanf(f, "%d\n", &maxval);

  int spacing = size / iwidth;
  int i, j;
  int val = 0;
  for (i = 0; i < iheight; i++)
  {
    for (j = 0; j < iwidth; j++)
    {
      fread(&val, 1, 1, f);
      fread(&val, 1, 1, f);
      fread(&val, 1, 1, f); // Nur B-Channel wird gelesen
      for (int il = 0; il < spacing; il++)
      {
        for (int jl = 0; jl < spacing; jl++)
        {
          int ib = i * spacing + il - imin + g;
          int jb = j * spacing + jl - jmin + g;
          if (jb >= g && jb < bwidth - g && ib >= g && ib < bheight - g)
            t[ib * bwidth + jb] = (double)val / 255 * 25.0; // Temperatur wird auf maximal 25.0 normiert
        }
      }
    }
  }
}
