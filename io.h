/**
 * @file io.h
 * @author Wigg, Seidler
 * @date 29 December 2019
 * @brief Methoden für Input und Output.
 */

#ifndef IO_H_
#define IO_H_

/**
 * @brief Liest die Parameter für die Simulation aus einer Datei ein.
 *
 * @param size Höhe/Breite des Gitters.
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
                   double *a, int *scenario, int *ostep, const char *filename);

/**
 * @brief Parallele Ausgabe eines Gitters als binärkodierte Portable Pixmap (PPM
 * P6).
 *
 * Schreibt des Gitterblock des Prozesses parallel in eine Outputdatei unter
 * Verwendung von MPI_File_write_at. Falls der Rang des Prozesses 0 ist, wird
 * zudem der Header geschrieben.
 *
 * @param t Gitter.
 * @param size Höhe/Breite des Gitters.
 * @param bwidth Breite des Gitterblocks.
 * @param bheight Höhe des Gitterblocks.
 * @param imin Index i der oberen linken Zelle des Gitterblocks.
 * @param jmin Index j der oberen linken Zelle des Gitterblocks.
 * @param g Breite der Geisterzonen.
 * @param filename Name der Output-Datei.
 * @param rank Rang des Prozesses, auf dem der Gitterblock liegt.
 */

void printPPM(double *t, int size, int bwidth, int bheight, int imin, int jmin,
              int g, const char *filename, int rank);

#endif  // IO_H_