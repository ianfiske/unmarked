
#ifndef _unmarked_DISTR_H
#define _unmarked_DISTR_H

#include <Rmath.h>
#include "utils.h"

// Zero-inflated Poisson
double dzip(int x, double lambda, double psi);

double N_density(int mixture, int k, double lambda, double log_alpha);

#endif
