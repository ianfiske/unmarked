#include "nll_pcountOpen.h"
#include "tranprobs.h"
#include "distr.h"


using namespace Rcpp ;



SEXP nll_pcountOpen( SEXP y_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xp_, SEXP Xiota_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_p_, SEXP beta_iota_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_) {
  int lk = as<int>(lk_);
  Rcpp::IntegerVector N = seq_len(lk)-1;
  int M = as<int>(M_);
  int J = as<int>(J_);
  int T = as<int>(T_);
  arma::imat ym = as<arma::imat>(y_);
  arma::mat Xlam = as<arma::mat>(Xlam_);
  arma::mat Xgam = as<arma::mat>(Xgam_);
  arma::mat Xom = as<arma::mat>(Xom_);
  arma::mat Xp = as<arma::mat>(Xp_);
  arma::mat Xiota = as<arma::mat>(Xiota_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_gam = as<arma::colvec>(beta_gam_);
  arma::colvec beta_om = as<arma::colvec>(beta_om_);
  arma::colvec beta_p = as<arma::colvec>(beta_p_);
  arma::colvec beta_iota = as<arma::colvec>(beta_iota_);
  double log_alpha = as<double>(log_alpha_);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offset_);
  arma::colvec Xgam_offset = as<arma::colvec>(Xgam_offset_);
  arma::colvec Xom_offset = as<arma::colvec>(Xom_offset_);
  arma::colvec Xp_offset = as<arma::colvec>(Xp_offset_);
  arma::colvec Xiota_offset = as<arma::colvec>(Xiota_offset_);
  std::string mixture = as<std::string>(mixture_);
  std::string dynamics = as<std::string>(dynamics_);
  std::string fix = as<std::string>(fix_);
  std::string go_dims = as<std::string>(go_dims_);
  bool immigration = as<bool>(immigration_);
  arma::imat I = as<arma::imat>(I_);
  arma::imat I1 = as<arma::imat>(I1_);
  Rcpp::List Ib(Ib_);
  Rcpp::List Ip(Ip_);
  int nrI = I.n_rows;
  int nrI1 = I1.n_rows;
  double alpha=0.0, psi=0.0;
  if(mixture=="NB")
    alpha = exp(log_alpha);
  else if(mixture=="ZIP")
    psi = 1.0/(1.0+exp(-log_alpha));
  Rcpp::IntegerVector first(first_);
  Rcpp::IntegerVector last(last_);
  arma::imat ytna = as<arma::imat>(ytna_); // y[i,,t] are all NA
  arma::imat ynam = as<arma::imat>(yna_);  // y[i,j,t] is NA
  arma::imat delta = as<arma::imat>(delta_);
  // linear predictors
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::mat omv = arma::ones<arma::colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend")) {
    if((dynamics == "ricker")  || (dynamics == "gompertz"))
        omv = exp(Xom*beta_om + Xom_offset);
    else if((dynamics == "constant")  || (dynamics == "autoreg") || (dynamics == "notrend"))
        omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
  }
  omv.reshape(T-1, M);
  arma::mat om = arma::trans(omv);
  arma::mat gam = arma::zeros<arma::mat>(M, T-1);
  if(dynamics == "notrend") {
    arma::mat lamMat = arma::repmat(lam, 1, T-1);
    gam = (1-om) % lamMat;
  } else {
    if(fix != "gamma") {
      arma::mat gamv = exp(Xgam*beta_gam + Xgam_offset);
      gamv.reshape(T-1, M);
      gam = arma::trans(gamv);
    }
  }
  arma::mat pv = 1.0/(1.0+exp(-1*(Xp*beta_p + Xp_offset)));
  pv.reshape(J*T, M);
  arma::mat pm = trans(pv);
  //Immigration
  arma::mat iotav = arma::zeros<arma::colvec>(M*(T-1));
  if(immigration) {
    iotav = exp(Xiota*beta_iota + Xiota_offset);
  }
  iotav.reshape(T-1, M);
  arma::mat iota = arma::trans(iotav);
  // format matrices as cubes (shouldn't be done in likelihood)
  arma::icube y(M,J,T);
  arma::cube p(M,J,T);
  arma::icube yna(M,J,T);
  for(int q=0; q<(M*J*T); q++) {
    y(q) = ym(q);
    p(q) = pm(q);
    yna(q) = ynam(q);
  }
  // initialize
  double ll=0.0;
  double ll_i=0.0;
  int first_i=0;
  int last_i=0;
  int first1=0;
  arma::colvec g_star = arma::ones<arma::colvec>(lk);
  arma::mat g3 = arma::zeros<arma::mat>(lk, lk);
  arma::mat g3_d = arma::zeros<arma::mat>(lk, lk);
  arma::colvec g1_t = arma::zeros<arma::colvec>(lk);
  arma::colvec g1_t_star = arma::zeros<arma::colvec>(lk);
  arma::colvec g1 = arma::zeros<arma::colvec>(lk);
  arma::colvec g2 = arma::zeros<arma::colvec>(lk);
  arma::colvec g1_star = arma::zeros<arma::colvec>(lk);
  arma::cube g3_t = arma::zeros<arma::cube>(lk,lk,T-1);

  // shouldn't be done in likelihood
  for(int i=0; i<M; i++) {
    if(first[i]==1) {
      first1=i; // site with no missing values at t=1
      break;
    }
  }

  // compute g3 if there are no covariates of omega/gamma
  if(go_dims == "scalar") {
    if(dynamics=="constant" || dynamics=="notrend")
      tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,0), om(first1,0));
    else if(dynamics=="autoreg")
      tp2(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="trend")
      tp3(g3, lk, gam(first1,0), iota(first1,0));
    else if(dynamics=="ricker")
      tp4(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="gompertz")
      tp5(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
  } else if(go_dims == "rowvec") {
    for(int t=0; t<(T-1); t++) {
      if(ytna(first1,t)==1) { // FIXME: this is not generic!
	continue;
      }
      if(dynamics=="constant" || dynamics=="notrend") {
	tp1(g3_t.slice(t), nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,t), om(first1,t));
      }
      else if(dynamics=="autoreg") {
	tp2(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
    }
      else if(dynamics=="trend")
	tp3(g3_t.slice(t), lk, gam(first1,t), iota(first1,t));
      else if(dynamics=="ricker")
	tp4(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
      else if(dynamics=="gompertz")
	tp5(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
    }
  }

  // loop over sites
  for(int i=0; i<M; i++) {
    first_i = first[i]-1; // remember 0=1st location in C++
    last_i = last[i]-1;
    g_star.ones();
    if(last_i > first_i) {
      // loop over time periods in reverse order, up to second occasion
      for(int t=last_i; t>first_i; t--) {
	if(ytna(i,t)==1) {
	  continue; //
	}
	g1_t.zeros();
	// loop over possible value of N at time t
	for(int k=0; k<lk; k++) {
	  for(int j=0; j<J; j++) {
	    if(yna(i,j,t)==0) {
	      g1_t(k) += Rf_dbinom(y(i,j,t), k, p(i,j,t), true);
	    }
	  }
	  g1_t(k) = exp(g1_t(k));
	  g1_t_star(k) = g1_t(k) * g_star(k);
	}
	// computes transition probs for g3
	if(go_dims == "matrix") {
	  g3.zeros();
	  if(dynamics=="constant" || dynamics=="notrend") {
	    //	    tp1(g3, lk, gam(i,t-1), om(i,t-1));
	    tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(i,t-1), om(i,t-1));
	  }
	  else if(dynamics=="autoreg") {
	    tp2(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	  }
	  else if(dynamics=="trend")
	    tp3(g3, lk, gam(i,t-1), iota(i,t-1));
	  else if(dynamics=="ricker")
	    tp4(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	  else if(dynamics=="gompertz")
	    tp5(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	} else if(go_dims == "rowvec") {
	  g3 = g3_t.slice(t-1);
	}
	int delta_it = delta(i,t);
	// matrix multiply transition probs over time gaps
	if(delta_it>1) {
	  g3_d = g3;
	  for(int d=1; d<delta_it; d++) {
	    g3_d = g3_d * g3;
	    }
	  g_star = g3_d * g1_t_star;
	} else
	  g_star = g3 * g1_t_star;
      }
    }
    ll_i=0.0;
    int delta_i0 = delta(i,0);
    g1.zeros();
    for(int k=0; k<lk; k++) { // loop over possible values of N
      for(int j=0; j<J; j++) {
	if(yna(i,j,first_i)==0) {
	  g1(k) += Rf_dbinom(y(i,j,first_i), k, p(i,j,first_i), true);
	}
      }
      g1(k) = exp(g1(k));
      if(delta_i0>1)
	g1_star(k) = g1(k) * g_star(k);
      if(mixture=="P")
	g2(k) = Rf_dpois(k, lam(i), false);
      else if(mixture=="NB")
        g2(k) = dnbinom_mu(k, alpha, lam(i), false);
      else if(mixture=="ZIP")
	g2(k) = dzip(k, lam(i), psi);
      if(delta_i0==1)
	ll_i += g1(k) * g2(k) * g_star(k);
    }
    if(delta_i0>1) {
      g3_d = g3;
      for(int d=0; d<delta_i0; d++) {
	g3_d = g3_d * g3;
      }
      g_star = g3_d * g1_star;
      ll_i = arma::dot(g2, g_star);
    }
    ll += log(ll_i + DBL_MIN);
  }
  return wrap(-ll);
}
