/**
 * @file io.h
 * @author Wigg, Seidler
 * @date 29 December 2019
 * @brief Methoden für Input und Output.
 */

#ifndef IO_H_
#define IO_H_

/**
 * @brief Ausgabe des Gitters als Portable Pix Map (PPM).
 *
 * Schreibt das Gitter t als PPM in die Datei filename.
 *
 * @param t Gitter.
 * @param size Seitenlänge des Gitters.
 * @param filename Dateiname der PPM.
 */
void printResult(double *t, int size, char *filename);

#endif  // IO_H_