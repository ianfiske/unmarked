#include "detfuns.h"
#include <Rmath.h>

// Half-normal detection function for point-transects
void grhn(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double sigma = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i]*x[i] / (2*sigma*sigma)) * x[i];
  }
}



