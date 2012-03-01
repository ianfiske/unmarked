#include "nll_distsamp.h"
#include "detfuns.h"

using namespace Rcpp ;

SEXP nll_distsamp( SEXP y_, SEXP X_, SEXP V_, SEXP beta_lam_, SEXP beta_sig_, SEXP X_offset_, SEXP V_offset_, SEXP A_, SEXP a_, SEXP u_, SEXP output_, SEXP db_ ) {

  arma::imat y = as<arma::imat>(y_);
  arma::mat X = as<arma::mat>(X_);
  arma::mat V = as<arma::mat>(V_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_sig = as<arma::colvec>(beta_sig_);
  arma::colvec X_offset = as<arma::colvec>(X_offset_);
  arma::colvec V_offset = as<arma::colvec>(V_offset_);
  arma::colvec A = as<arma::colvec>(A_);
  arma::mat a = as<arma::mat>(a_);
  arma::mat u = as<arma::mat>(u_);
  std::string output = as<std::string>(string_);
  Rcpp::NumericVector db(db_);

  int R = y.n_rows;
  int J = y.n_cols;
  arma::colvec lam = exp(X*beta_lam + X_offset);
  arma::colvec sig = exp(V*beta_sig + V_offset);
  // need to add output option to model density instead of abundance
  if(output=="density")
    lam = lam % A; // double check %

  // Integration settings given to Rdqags
  double epsabs = 0.01; // should be specific to measurement units
  double epsrel = 0.01; // should be specific to measurement units
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
      ex = sig(i);
      a = db[j];
      b = db[j+1];
      double result = 0.0;
      double abserr = 0.0;
      int neval = 0;
      int ier = 0;
      Rdqags(grhn, ex, &a, &b, &epsabs, &epsrel, &result, &abserr,
	     &neval, &ier, &limit, &lenw, &last, &iwork, &work);
      cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
      mu = lam(i)*cp + DOUBLE_XMIN;
      ll += Rf_dpois(y(i,j), lam(i)*cp, true);
    }
  }
  return wrap(-ll);
}
