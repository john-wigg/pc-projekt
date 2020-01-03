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
 * @param size Seitenlänge des Gitters.
 * @param filename Dateiname der PPM.
 */
void printPPMP3(double *t, int size, char *filename);

/**
 * @brief Ausgabe des Gitters als binärkodierte Portable Pix Map (PPM).
 *
 * Schreibt das Gitter t als binärkodierte PPM in die Datei filename.
 *
 * @param t Gitter.
 * @param size Seitenlänge des Gitters.
 * @param filename Dateiname der PPM.
 */
void printPPMP6(double *t, int size, char *filename);

void readGraymap(int **idi, int **idj, double **t, int size, int bheight,
                 int bwidth, int imin, int jmin, int g, char *filename);

#endif  // IO_H_