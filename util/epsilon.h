#ifndef EPSILON_H
#define EPSILON_H

/* Utility functions for epsilon equality testing of doubles. Two doubles
 * are considered equal if they are within EPSILON defined in epsilon.c.
 * For epsilon_less_or_equal use !epsilon_greater.
 * For epsilon_greater_or_equal use !epsilon_less.
 */
int epsilon_equal(double a, double b);   // a == b.
int epsilon_less(double a, double b);    // a <  b.
int epsilon_greater(double a, double b); // a >  b.

#endif // EPSILON_H
