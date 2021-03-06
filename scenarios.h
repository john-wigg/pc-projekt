/**
 * @file scenarios.h
 * @author Wigg
 * @date 30 January 2020
 * @brief Funktionen zur Initialisierung der Gitterblöcke mit verschiendenen
 * Anfangs-/ Randbedingungen.
 */

#ifndef SCENARIOS_H_
#define SCENARIOS_H_

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

/**
 * @brief Initialisiert einen Gitterblock mit einem heißen "Fleck" in der Mitte
 * des Feldes.
 *
 * Siehe initTopLeftHotBotRightCold() für Parameterliste.
 */
void initHotSpot(double *t, int size, int bheight, int bwidth, int imin,
                 int jmin, int g);

/**
 * @brief Initialisiert einen Gitterblock mit 25.0 auf dem rechten Rand.
 *
 * Siehe initTopLeftHotBotRightCold() für Parameterliste.
 */
void initRightHot(double *t, int size, int bheight, int bwidth, int imin,
                  int jmin, int g);

#endif  // SCENARIOS_H_