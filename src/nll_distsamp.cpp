#include "nll_distsamp.h"
#include "detfuns.h"

SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP a_, SEXP u_, SEXP db_ ) {

  Rcpp::IntegerMatrix y(y_);
  Rcpp::NumericVector lam(lam_);
  Rcpp::NumericVector sig(sig_);
  Rcpp::NumericMatrix a(a_);
  Rcpp::NumericMatrix u(u_);
  Rcpp::NumericVector db(db_);

  int R = y.nrow();   //y.n_rows;
  int J = y.ncol();   // y.n_cols;

  // Integration settings given to Rdqags
  // TODO: Allow for user-defined settings
  double epsabs = 0.01; // should be specific to measurement units
  double epsrel = 0.00012; // should be specific to measurement units
  int limit = 100;
  int lenw = 400;
  int last = 0;
  int iwork = 100;
  double work = 400.0;

  double mu = 0.0;
  double ll = 0.0;

  for(int i=0; i<R; i++) {
    for(int j=0; j<J; j++) {
      double cp = 0.0;
      void *ex;
      ex = &sig(i);
      double lower = db[j];
      double upper = db[j+1];
      double result = 0.0;
      double abserr = 0.0;
      int neval = 0;
      int ier = 0;
      Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result, &abserr,
	     &neval, &ier, &limit, &lenw, &last, &iwork, &work);
      /* add error checking/handling here */
      cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
      mu = lam(i)*cp + DOUBLE_XMIN;
      ll += Rf_dpois(y(i,j), lam(i)*cp, true);
    }
  }
  return Rcpp::wrap(-ll);
}
