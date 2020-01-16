/**
 * @file io.h
 * @author Wigg, Seidler
 * @date 29 December 2019
 * @brief Methoden für Input und Output.
 */

#ifndef IO_H_
#define IO_H_

/**
 * @brief Ausgabe des Gitters als ASCII-kodierte Portable Pix Map (PPM).
 *
 * Schreibt das Gitter t als ASCII-kodierte PPM in die Datei filename.
 *
 * @param t Gitter.
 * @param size Anzahl Gitterpunkte in eine Richtung.
 * @param filename Dateiname der PPM.
 */
void printPPMP3(double *t, int size, const char *filename);

/**
 * @brief Ausgabe des Gitters als binärkodierte Portable Pix Map (PPM).
 *
 * Schreibt das Gitter t als binärkodierte PPM in die Datei filename.
 *
 * @param t Gitter.
 * @param size Anzahl Gitterpunkte in eine Richtung.
 * @param filename Dateiname der PPM.
 */
void printPPMP6(double *t, int size, const char *filename);

/**
 * @brief Liest die Parameter für die Simulation aus einer Datei ein.
 *
 * @param size Anzahl Gitterpunkte in eine Richtung.
 * @param iter Anzahl Iterationen.
 * @param g Breite der Überschneidungsbereiche.
 * @param adj Sicherheitsfaktor des Zeitschrittes (<1.0 für stabiles
 * Euler-Verfahren)
 * @param alpha Thermische Diffusivität.
 * @param a Seitenlänge des Gitters.
 * @param scenario Simulationsszenario. Mögliche Werte: TOPLEFTHOTBOTRIGHTCOLD,
 * RIGHTHOT, HOTSPOT, oder Dateiname um eine PPM-Datei einzulesen.
 * @param filename Dateiname der Inputdatei.
 */
void readInputFile(int *size, int *iter, int *g, double *adj, double *alpha,
                   double *a, int *scenario, const char *filename);

void printPPMP6MPI(double *t, int size, int bwidth, int bheight, int imin,
                   int jmin, int g, const char *filename, int rank);

#endif  // IO_H_