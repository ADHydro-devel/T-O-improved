#include <math.h>
#include "epsilon.h"

// Epsilon is the larger of 10^-10 or ten orders of magnitude smaller than the number being checked because large doubles might have less resolution than ten
// digits past the decimal point, but they will always have ten significant digits.
#define EPSILON(x) ((1.0 >= fabs(x)) ? (1.0e-10) : (fabs(x) * 1.0e-10))

/* Comment in .h file. */
int epsilon_equal(double a, double b)
{
  return !epsilon_less(a, b) && !epsilon_greater(a, b);
}

/* Comment in .h file. */
int epsilon_less(double a, double b)
{
  return a < b - EPSILON(b);
}

/* Comment in .h file. */
int epsilon_greater(double a, double b)
{
  return a > b + EPSILON(b);
}
