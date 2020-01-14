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
 * ihm, so wird der Tag \p T_SEND_E verwendet.
 */
enum {
  T_SEND_N,
  T_SEND_S,
  T_SEND_E,
  T_SEND_W,
  T_SEND_NW,
  T_SEND_NE,
  T_SEND_SW,
  T_SEND_SE
};

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