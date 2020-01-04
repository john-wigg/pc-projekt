#ifndef SCENARIOS_H_
#define SCENARIOS_H_

#include "scenarios.h"

/**
 * @brief Initialisiert einen Gitterblock mit Randbedingungen.
 *
 * Initialisiert den Gitterblock mit 0.0 für alle inneren Punkte, mit 25.0 am
 * oberen und am linken Rand und mit -25.0 am unteren und am rechten Rand des
 * Gesamtgitters.
 *
 * @param t Gitterblock.
 * @param size Seitenlänge des Gesamtgitters.
 * @param bsize Seitenlänge der Gitterblocks.
 * @param imin Index i des oberen linken Punktes des Blocks.
 * @param jmin Index j des oberen linken Punktes des Blocks.
 * @param g Überschneidungsbreite der Blöcke.
 */
void initTopLeftHotBotRightCold(double *t, int size, int bheight, int bwidth,
                                int imin, int jmin, int g);

void initHotSpot(double *t, int size, int bheight, int bwidth, int imin,
                 int jmin, int g);

void initRightHot(double *t, int size, int bheight, int bwidth, int imin,
                  int jmin, int g);

void initFromImage(double *t, int size, int bheight, int bwidth, int imin,
                   int jmin, int g, char *filename);

#endif  // SCENARIOS_H_