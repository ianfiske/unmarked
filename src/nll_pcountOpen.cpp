#include "nll_pcountOpen.h"

using namespace Rcpp ;

SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP ytr_, SEXP yr_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_ ) {
  int lk = as<int>(lk_);
  int M = as<int>(M_);
  int J = as<int>(J_);
  int T = as<int>(T_);
  arma::imat ym = as<arma::imat>(y_);
  arma::mat Xlam = as<arma::mat>(Xlam_);
  arma::mat Xgam = as<arma::mat>(Xgam_);
  arma::mat Xom = as<arma::mat>(Xom_);
  arma::mat Xp = as<arma::mat>(Xp_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_gam = as<arma::colvec>(beta_gam_);
  arma::colvec beta_om = as<arma::colvec>(beta_om_);
  arma::colvec beta_p = as<arma::colvec>(beta_p_);
  double log_alpha = as<double>(log_alpha_);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offset_);
  arma::colvec Xgam_offset = as<arma::colvec>(Xgam_offset_);
  arma::colvec Xom_offset = as<arma::colvec>(Xom_offset_);
  arma::colvec Xp_offset = as<arma::colvec>(Xp_offset_);
  double alpha = exp(log_alpha);
  std::string mixture = as<std::string>(mixture_);
  Rcpp::NumericVector first(first_);
  Rcpp::NumericVector last(last_);
  arma::imat ytr = as<arma::imat>(ytr_); // y[i,,t] are not all NA
  arma::imat yrm = as<arma::imat>(yr_);  // y[i,j,t] is not NA
  arma::imat delta = as<arma::imat>(delta_);
  // linear predictors
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec gamv = exp(Xgam*beta_gam + Xgam_offset);
  arma::colvec omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
  arma::colvec pv = 1.0/(1.0+exp(-1*(Xp*beta_p + Xp_offset)));
  // Put vectors in row-major matrices
  gamv.reshape(T-1, M);
  arma::mat gam = trans(gamv);
  omv.reshape(T-1, M);
  arma::mat om = trans(omv);
  pv.reshape(T, M);
  arma::mat pm = trans(pv);
  // format matrices as cubes
  arma::cube y(M,J,T);
  arma::cube p(M,J,T);
  arma::icube yr(M,J,T);
  for(int q=0; q<M*J*T; q++) {
    y(q) = ym(q);
    p(q) = pm(q);
    yr(q) = yrm(q);
  }
  // initialize
  double ll=0.0;
  double ll_i=0.0, g1=0.0, g2=0.0;
  int first_i=0, last_i=0;
  arma::colvec g_star = arma::ones<arma::colvec>(lk);
  arma::cube g3 = arma::zeros<arma::cube>(lk, lk, T-1);
  arma::colvec g1_t = arma::zeros<arma::colvec>(lk);
  arma::colvec g1_t_star = arma::zeros<arma::colvec>(lk);
  // loop over sites
  for(int i=0; i<M; i++) {
    first_i = first[i]-1; // remember 0=1st location in C++
    last_i = last[i]-1;
    g_star.ones();
    g3.zeros();
    if(last_i > first_i) {
      // loop over time periods in reverse order, up to second occasion
      for(int t=last_i; t>first_i; t--) {
	if(ytr(i,t)==0) {
	  continue; // FIXME: g3 needs to handle time gap
	}
	g1_t.zeros();
	// loop over possible value of N at time t
	for(int n2=0; n2<lk; n2++) {
	  if(J==1)
	    g1_t(n2) = Rf_dbinom(y(i,0,t), n2, p(i,0,t), false);
	  else {
	    for(int j=0; j<J; j++) {
	      if(yr(i,j,t)==1) {
		g1_t(n2) += Rf_dbinom(y(i,j,t), n2, p(i,j,t), true);
		g1_t(n2) = exp(g1_t(n2));
	      }
	    }
	  }
	  g1_t_star(n2) = g1_t(n2) * g_star(n2);
	}
	// computes transition probs for g3.slice(t-1)
	tp(g3, lk, gam(i,t-1), om(i,t-1), t-1);
	g_star = g3.slice(t-1) * g1_t_star;
      }
    }
    g1=0.0;
    g2=0.0;
    ll_i=0.0;
    for(int n1=0; n1<lk; n1++) { // loop over possible values of N
      if(J==1)
	g1 = Rf_dbinom(y(i,0,first_i), n1, p(i,0,first_i), false);
      else {
	for(int j=0; j<J; j++) {
	  if(yr(i,j,first_i)==1) {
	    g1 += Rf_dbinom(y(i,j,first_i), n1, p(i,j,first_i), true);
	    g1 = exp(g1);
	  }
	}
      }
      if(mixture=="P")
	g2 = Rf_dpois(n1, lam(i), false);
      else
        g2 = dnbinom_mu(n1, alpha, lam(i), false);
      //      if(delta(i,0)==1)
	ll_i += g1 * g2 * g_star(n1);
	//      else if(delta(i,0)>1) {
	// FIXME: need to recompute g3 here, call it g3_0
	// g1_t_star(n1) = g1 * g_star(n1);
	// g_star(n1) = trans(g3_0) * g1_t_star; // this is wrong
	// ll_i += g2 * g_star(n1);
	//}
    }
    ll += log(ll_i);
  }
  return wrap(-ll);
}
