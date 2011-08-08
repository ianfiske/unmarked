#include "nll_pcount.h"

using namespace Rcpp ;

SEXP nll_pcount( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_lamR, SEXP beta_pR, SEXP log_alphaR, SEXP X_offsetR, SEXP V_offsetR, SEXP naMatR, SEXP lkR, SEXP mixtureR ) {
  arma::mat y = as<arma::mat>(yR);
  arma::mat X = as<arma::mat>(Xr);
  arma::mat V = as<arma::mat>(Vr);
  arma::colvec beta_lam = as<arma::colvec>(beta_lamR);
  arma::colvec beta_p = as<arma::colvec>(beta_pR);
  double log_alpha = as<double>(log_alphaR);
  double alpha = exp(log_alpha);
  arma::colvec X_offset = as<arma::colvec>(X_offsetR);
  arma::colvec V_offset = as<arma::colvec>(V_offsetR);
  Rcpp::LogicalMatrix naMat(naMatR);
  std::string mixture = as<std::string>(mixtureR);
  int lk = as<int>(lkR);
  int R = X.n_rows;
  int J = y.n_cols;
  arma::colvec lam = exp(X*beta_lam + X_offset);
  Rcpp::NumericMatrix pmat(R,J);
  arma::colvec logit_p = V*beta_p + V_offset;
  arma::colvec p = 1.0/(1.0+exp(-logit_p));
  double ll=0.0;
  Rcpp::NumericVector f(lk);
  Rcpp::NumericVector g(lk);
  int z=0, zi=0;
  for(int i=0; i<R; i++) {
    zi=i*J;
    for(int k=0; k<lk; k++) {
      if(mixture=="P")
	f(k) = Rf_dpois(k, lam(i), false);
      else
        f(k) = dnbinom_mu(k, alpha, lam(i), false);
      g(k) = 1.0;
      for(int j=0; j<J; j++) {
	// if(k >= y(i,j))
	z = zi + j;
	if(!naMat(i,j))
	  g(k) *= Rf_dbinom(y(i,j), k, p(z), false);
	// else
	// g(k) = 0.0;
      }
    }
    ll += log(sum(f*g));
  }
  return wrap(-ll);

}
