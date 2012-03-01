#include "detfuns.h"

// Half-normal detection function for point-transects
void grhn(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double sigma = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i]*x[i] / (2*sigma*sigma)) * x[i];
  }
}




/*
  void f(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex;
  double sigma = v[0];
  for(int i=0; i<n; i++) {
  x[i] = exp(-(x[i]*x[i]) / (2*sigma*sigma));
  }
  }



    double sigma = as<double>(sigma_);
    double a = as<double>(a_);
    double b = as<double>(b_);
    void *ex;
    ex = &sigma;
    double epsabs = 0.001;
    double epsrel = 0.001;
    double result = 0.0;
    double abserr = 0.0;
    int neval = 0;
    int ier = 0;
    int limit = 100;
    int lenw = 400;
    int last = 0;
    int iwork=100;
    double work=400.0;
    Rdqags(f, ex, &a, &b, &epsabs, &epsrel, &result, &abserr, &neval, &ier,
           &limit, &lenw, &last, &iwork, &work);
    return Rcpp::List::create(result, abserr, last, ier);

*/
