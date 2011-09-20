#include "nll_pcountOpen.h"

using namespace Rcpp ;

SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP log_alpha_, SEXP naMat_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_ ) {
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
  Rcpp::LogicalMatrix naMat(naMat_);
  std::string mixture = as<std::string>(mixture_);
  Rcpp::NumericVector first(first_);
  Rcpp::NumericVector last(last_);
  int lk = as<int>(lk_);
  int R = y.n_rows;
  int T = y.n_cols;
  // offsets ignored
  arma::colvec lam = exp(Xlam*beta_lam);
  arma::colvec gamv = exp(Xgam*beta_gam);
  arma::colvec omv = 1.0/(1.0+exp(-1*(Xom*beta_om)));
  arma::colvec pv = 1.0/(1.0+exp(-1*(Xp*beta_p)));
  arma::mat gam = reshape(gamv, R, T-1, 1); // row-wise
  arma::mat om = reshape(omv, R, T-1, 1);
  arma::mat p = reshape(pv, R, T, 1);
  double ll=0.0;
  double ll_i, g1, g2;
  int Nmin=0;
  arma::colvec g_star(lk);
  arma::cube g3 = arma::cube(lk, lk, T-1);
  arma::colvec g1_t(lk);
  arma::colvec g1_t_star(lk);
  for(int i=0; i<R; i++) {
    int first_i = first[i]-1;
    int last_i = last[i]-1;
    ll_i=0.0;
    g_star.ones();
    for(int t=last_i; t>0; t--) { // last through 2nd occassion
      for(int n2=0; n2<lk; n2++) {
        g1_t(n2) = Rf_dbinom(y(i, t), n2, p(i, t), false);
        g1_t_star(n2) = g1_t(n2) * g_star(n2);
        for(int n1=0; n1<lk; n1++) {
          Nmin = std::min(n1, n2);
          for(int c=0; c<=Nmin; c++) {
	    g3(n2, n1, t-1) += exp(Rf_dbinom(c, n1, om(i,t-1), true) +
				   Rf_dpois(n2-c, gam(i,t-1), true));
	  }
	}
      }
      g_star = g3.slice(t-1) * g1_t_star;
    }
    g1=0.0;
    g2=0.0;
    for(int k=0; k<lk; k++) {
      g1 = Rf_dbinom(y(i, first_i), k, p(i, first_i), false);
      if(mixture=="P")
	g2 = Rf_dpois(k, lam(i), false);
      else
        g2 = dnbinom_mu(k, alpha, lam(i), false);
      ll_i += g1 * g2 * g_star(k);
    }
    ll += log(ll_i);
  }
  return wrap(-ll);
}
