#include "nll_pcountOpen.h"

using namespace Rcpp ;

SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP log_alpha_, SEXP X_offset_, SEXP V_offsetR, SEXP naMat_, SEXP lk_, SEXP mixture_ ) {
  arma::mat y = as<arma::mat>(y_);
  arma::mat Xlam = as<arma::mat>(Xlam_);
  arma::mat Xgam = as<arma::mat>(Xgam_);
  arma::mat Xom = as<arma::mat>(Xom_);
  arma::mat Xp = as<arma::mat>(Xp_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_gam = as<arma::colvec>(beta_gam_);
  arma::colvec beta_om = as<arma::colvec>(beta_om_);
  arma::colvec beta_p = as<arma::colvec>(beta_p_);
  double log_alpha = as<double>(log_alpha_);
  double alpha = exp(log_alpha);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offsetR_);
  arma::colvec Xgam_offset = as<arma::colvec>(Xgam_offsetR_);
  arma::colvec Xom_offset = as<arma::colvec>(Xom_offsetR_);
  arma::colvec Xp_offset = as<arma::colvec>(Xp_offsetR_);
  Rcpp::LogicalMatrix naMat(naMat_);
  std::string mixture = as<std::string>(mixture_);
  int lk = as<int>(lk_);
  int R = X.n_rows;
  int T = y.n_cols;
  arma::colvec lamv = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec gamv = exp(Xgam*beta_gam + Xgam_offset);
  arma::colvec omv = 1.0/(1.0+exp(-(Xom*beta_om + Xom_offset)));
  arma::colvec pv = 1.0/(1.0+exp(-(Xp*beta_p + Xp_offset)));
  arma::mat gam = reshape(gamv, R, T, 1);
  arma::mat om = reshape(omv, R, T, 1);
  arma::mat p = reshape(pv, R, T, 1);
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
