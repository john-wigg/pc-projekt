/**
 * @file blocks.h
 * @author Wigg
 * @date 29 December 2019
 * @brief Methoden zur Arbeit mit Gitterblöcken.
 */

#ifndef BLOCKS_H_
#define BLOCKS_H_

/**
 * @brief Tags zur  Zuordnung der Überschneidungsbereiche.
 *
 * Gibt an, in welche Richtung der Überschneidungsbereich gesendet wird. Sendet
 * zum Beispiel ein Block ihren Überschneidungsbereich an den Block links von
 * ihm, so wird der Tag \p T_SENDLEFT verwendet.
 */
enum { T_SENDUP, T_SENDDOWN, T_SENDLEFT, T_SENDRIGHT };

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
void initBlock(double *t, int size, int bheight, int bwidth, int imin, int jmin,
               int g);

/**
 * @brief Vertauscht die Pointer zweier Gitterblöcke.
 *
 * @param t1 Pointer zum ersten Gitterblock.
 * @param t2 Pointer zum zweiten Gitterblokc.
 */
void swap(double **t1, double **t2);

/**
 * @brief Aktualisiert einen Gitterpunkt.
 *
 * Führt einen einzelnen Euler-Vorwärts-Schritt für einen Punkt auf dem Gitter
 * aus und speichert das Ergebnis in ein neues Gitter.
 *
 * @param t1 Neues Gitter.
 * @param t2 Altes Gitter.
 * @param i Index i des Punktes auf dem Gitter.
 * @param j Index j des Punktes auf dem Gitter.
 * @param adjstep Angepasste Schrittweite des Euler-Schritts: \f$\alpha
 * \frac{\Delta t}{h^2}\f$.
 */
void updatePoint(double *t1, double *t2, int i, int j, double adjstep,
                 int size);

#endif  // BLOCKS_H_