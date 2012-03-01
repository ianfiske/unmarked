
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "detfuns.h"

// Half-normal detection function for point-transects
double grhn(int x, double sigma) {
  double p=0.0;
  p = exp(-x*x / (2*sigma*sigma)) * x
  return p;
}

