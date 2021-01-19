
#include "distr.h"

// Zero-inflated Poisson
double dzip(int x, double lambda, double psi) {
  double den=0.0;
  if(x==0)
    den = psi + (1-psi)*exp(-lambda);
  else if(x>0)
    den = (1-psi)*Rf_dpois(x, lambda, false);
  return den;
}

double N_density(int mixture, int x, double lambda, double alpha){
  switch(mixture){
    case 1: return Rf_dpois(x, lambda, false);
    case 2: return Rf_dnbinom_mu(x, exp(alpha), lambda, false);
    case 3: return dzip(x, lambda, inv_logit(alpha));
  }
  return 0.0;
}
