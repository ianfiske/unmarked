#include "nll_pcountOpen.h"

using namespace Rcpp ;

SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP naMat_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_ ) {
  arma::mat ym = as<arma::mat>(y_);
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
  Rcpp::NumericVector first(first_);
  Rcpp::NumericVector last(last_);
  int lk = as<int>(lk_);
  int R = X.n_rows;
  int T = y.n_cols;
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec gamv = exp(Xgam*beta_gam + Xgam_offset);
  arma::colvec omv = 1.0/(1.0+exp(-(Xom*beta_om + Xom_offset)));
  arma::colvec pv = 1.0/(1.0+exp(-(Xp*beta_p + Xp_offset)));
  arma::mat gam = reshape(gamv, R, T, 1); // row-wise
  arma::mat om = reshape(omv, R, T, 1);
  arma::mat p = reshape(pv, R, T, 1);
  double ll=0.0;
  double g1=0.0, g2=0.0;
  int Nmin=0;
  for(int i=0; i<R; i++) {
    int first_i = first[i]-1;
    int last_i = last[i]-1;
    int ll_i=0.0;
    arma::colvec g_star = arma::ones<arma::colvec>(lk);
    arma::cube g3 = arma::cube(lk, lk, T-1);
    arma::colvec g1_t(lk);
    for(int t=last_i; t>0; t--) { // last through 2nd occassion
      for(int N2=0; N2<lk; N2++) {
        g1_t(N2) = Rf_binom(y(i, t), N2, p(i, t), true);
        for(int N1=0; N1<lk; N1++) {
          Nmin = std::min(N1, N2);
          for(int c=0; c<=Nmin; c++) {
	    g3(N1, N2, t-1) += exp(Rf_dbinom(c, N1, om(i,t-1), true) +
				   Rf_dpois(N2-c, gam(i,t-1), true));
	  }
	}
      }
      g_star = g3_t.slice(t-1) * g1_t % g_star(N2);
    }
    for(int k=0; k<lk; k++) {
      g1 = Rf_dbinom(y(i, first_i), k, p(i, first_i), false);
      if(mixture=="P")
	g2 = Rf_dpois(k, lam(i), false);
      else
        g2 = dnbinom_mu(k, alpha, lam(i), false);
      ll_i += g1 * g2 * g_star[k];
    }
    ll += log(ll_i);
  }
  return wrap(-ll);

}
