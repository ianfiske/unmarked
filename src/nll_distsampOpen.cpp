#include "nll_distsampOpen.h"
#include "distr.h"
#include "detfuns.h"

using namespace Rcpp ;


SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xsig_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_sig_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP fix_, SEXP go_dims_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP a_, SEXP u_, SEXP db_ ) {

  Rprintf("check 1\n");

  int lk = as<int>(lk_);
  Rcpp::IntegerVector N = seq_len(lk)-1;
  int M = as<int>(M_);
  int J = as<int>(J_);
  int T = as<int>(T_);
  arma::imat ym = as<arma::imat>(y_);
  arma::imat yt = as<arma::imat>(yt_);
  arma::mat Xlam = as<arma::mat>(Xlam_);
  arma::mat Xgam = as<arma::mat>(Xgam_);
  arma::mat Xom = as<arma::mat>(Xom_);
  arma::mat Xsig = as<arma::mat>(Xsig_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_gam = as<arma::colvec>(beta_gam_);
  arma::colvec beta_om = as<arma::colvec>(beta_om_);
  arma::colvec beta_sig = as<arma::colvec>(beta_sig_);
  double log_alpha = as<double>(log_alpha_);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offset_);
  arma::colvec Xgam_offset = as<arma::colvec>(Xgam_offset_);
  arma::colvec Xom_offset = as<arma::colvec>(Xom_offset_);
  arma::colvec Xsig_offset = as<arma::colvec>(Xsig_offset_);
  std::string mixture = as<std::string>(mixture_);
  std::string dynamics = as<std::string>(dynamics_);
  std::string fix = as<std::string>(fix_);
  std::string go_dims = as<std::string>(go_dims_);
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

  Rprintf("check 2\n");

  arma::mat a = as<arma::mat>(a_);
  arma::mat u = as<arma::mat>(u_);
  Rcpp::NumericVector db(db_);

  // linear predictors
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec omv = arma::ones<arma::colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend"))
    omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
  omv.reshape(T-1, M);
  arma::mat om = arma::trans(omv);
  arma::mat gam = arma::zeros<arma::mat>(M, T-1);
  if(dynamics == "notrend") {
    arma::mat lamMat = arma::repmat(lam, 1, T-1);
    gam = (1-om) % lamMat;
  } else {
    if(fix != "gamma") {
      arma::colvec gamv = exp(Xgam*beta_gam + Xgam_offset);
      gamv.reshape(T-1, M);
      gam = arma::trans(gamv);
    }
  }
  arma::colvec sigv = exp(Xsig*beta_sig + Xsig_offset);
  sigv.reshape(T, M);
  arma::mat sig = trans(sigv);
  // format matrices as cubes (shouldn't be done in likelihood)
  arma::icube y(M,J,T);
  arma::icube yna(M,J,T);
  for(int q=0; q<(M*J*T); q++) {
    y(q) = ym(q);
    yna(q) = ynam(q);
  }

  Rprintf("check 3\n");

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

  // Integration settings given to Rdqags
  double epsabs = 0.01; // should be specific to measurement units
  double epsrel = 0.01; // should be specific to measurement units
  int limit = 100;
  int lenw = 400;
  int last2 = 0;
  int iwork = 100;
  double work = 400.0;

  // shouldn't be done in likelihood
  for(int i=0; i<M; i++) {
    if(first[i]==1) {
      first1=i; // site with no missing values at t=1
      break;
    }
  }

  double cp = 0.0;
  double cpsum = 0.0;
  double part1 = 0.0;
  double part2 = 0.0;
  double part3 = 0.0;

  // compute g3 is there are no covariates of omega/gamma
  if(go_dims == "scalar") {
    if(dynamics=="constant" || dynamics=="notrend") {
      tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,0), om(first1,0));
    }
    else if(dynamics=="autoreg") {
      tp2(g3, lk, gam(first1,0), om(first1,0));
    }
    else if(dynamics=="trend")
      tp3(g3, lk, gam(first1,0));
  } else if(go_dims == "rowvec") {
    for(int t=0; t<(T-1); t++) {
      if(ytna(first1,t)==1) { // FIXME: this is not generic!
	continue;
      }
      if(dynamics=="constant" || dynamics=="notrend") {
	tp1(g3_t.slice(t), nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,t), om(first1,t));
      }
      else if(dynamics=="autoreg") {
	tp2(g3_t.slice(t), lk, gam(first1,t), om(first1,t));
      }
      else if(dynamics=="trend")
	tp3(g3_t.slice(t), lk, gam(first1,t));
    }
  }

  Rprintf("check 4\n");

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

//	Rprintf("check 5\n");

	for(int k=yt(i,t); k<lk; k++) { // note first index
	  for(int j=0; j<=J; j++) { // note <=J not <J
	    part1 = lgamma(k+1); // compute outside likelihood, or ignore
	    cp = 0.0;
	    cpsum = 0.0;
	    if(j<J) {
	      void *ex;
	      ex = &sig(i,t);
	      double lower = db[j];
	      double upper = db[j+1];
	      double result = 0.0;
	      double abserr = 0.0;
	      int neval = 0;
	      int ier = 0;
	      Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
		     &abserr, &neval, &ier, &limit, &lenw, &last2, &iwork,
		     &work);
	      /* add error checking/handling here */
	      cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	      cpsum = cpsum+cp;
	      part2 = lgamma(y(i,j,t)+1); // compute outside likelihood
	      part3 = log(cp + DOUBLE_XMIN) * y(i,j,t);
	    } else {
	      cp = 1 - cpsum;
	      part2 = lgamma(k-yt(i,t)+1);
	      part3 = log(cp + DOUBLE_XMIN) * (k - yt(i,t));
	    }
	    g1_t(k) += part1 - part2 + part3;
	  }
	  g1_t(k) = exp(g1_t(k));
	  g1_t_star(k) = g1_t(k) * g_star(k);
	}

//	Rprintf("check 6\n");

	// computes transition probs for g3
	if(go_dims == "matrix") {
	  g3.zeros();
	  if(dynamics=="constant" || dynamics=="notrend") {
	    //	    tp1(g3, lk, gam(i,t-1), om(i,t-1));
	    tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(i,t-1), om(i,t-1));
	  }
	  else if(dynamics=="autoreg") {
	    tp2(g3, lk, gam(i,t-1), om(i,t-1));
	  }
	  else if(dynamics=="trend")
	    tp3(g3, lk, gam(i,t-1));
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

//    Rprintf("check 7\n");

    ll_i=0.0;
    int delta_i0 = delta(i,0);
    g1.zeros();
    for(int k=yt(i,first_i); k<lk; k++) { // loop over possible values of N
      for(int j=0; j<=J; j++) {
	part1 = lgamma(k+1); // compute outside likelihood, or ignore
	cp = 0.0;
	cpsum = 0.0;
	if(j<J) {
	  void *ex;
	  ex = &sig(i,first_i);
	  double lower = db[j];
	  double upper = db[j+1];
	  double result = 0.0;
	  double abserr = 0.0;
	  int neval = 0;
	  int ier = 0;
	  Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
		 &abserr, &neval, &ier, &limit, &lenw, &last2, &iwork,
		 &work);
	  /* add error checking/handling here */
	  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	  cpsum = cpsum+cp;
	  part2 = lgamma(y(i,j,first_i)+1); // compute outside likelihood
	  part3 = log(cp + DOUBLE_XMIN) * y(i,j,first_i);
	} else {
	  cp = 1 - cpsum;
	  part2 = lgamma(k-yt(i,first_i)+1);
	  part3 = log(cp + DOUBLE_XMIN) * (k - yt(i,first_i));
	}
	g1(k) += part1 - part2 + part3;
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

    Rprintf("Check 8\n");
    Rprintf("site : i%\n", i);

    if(delta_i0>1) {
      g3_d = g3;
      for(int d=0; d<delta_i0; d++) {
	g3_d = g3_d * g3;
      }
      g_star = g3_d * g1_star;
      ll_i = arma::dot(g2, g_star);
    }
    ll += log(ll_i + DOUBLE_XMIN);
  }

  Rprintf("Made it\n");

  return wrap(-ll);
}
