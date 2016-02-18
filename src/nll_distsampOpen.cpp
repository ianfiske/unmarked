#include "nll_distsampOpen.h"
#include "tranProbs.h"
#include "distr.h"
#include "detfuns.h"

using namespace Rcpp;



SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, SEXP Xsig_, SEXP Xiota_, SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, SEXP beta_sig_, SEXP beta_iota_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP Xiota_offset_, SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, SEXP M_, SEXP J_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP survey_, SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP scale_, SEXP a_, SEXP u_, SEXP db_, SEXP nintervals_, SEXP keyfun_) {

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
  arma::mat Xiota = as<arma::mat>(Xiota_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_gam = as<arma::colvec>(beta_gam_);
  arma::colvec beta_om = as<arma::colvec>(beta_om_);
  arma::colvec beta_sig = as<arma::colvec>(beta_sig_);
  arma::colvec beta_iota = as<arma::colvec>(beta_iota_);
    double log_alpha = as<double>(log_alpha_);
//    double scale = Rcpp::as<double>(scale_);
//  double scale = as<double>(scale_);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offset_);
  arma::colvec Xgam_offset = as<arma::colvec>(Xgam_offset_);
  arma::colvec Xom_offset = as<arma::colvec>(Xom_offset_);
  arma::colvec Xsig_offset = as<arma::colvec>(Xsig_offset_);
  arma::colvec Xiota_offset = as<arma::colvec>(Xiota_offset_);
  std::string mixture = as<std::string>(mixture_);
  std::string dynamics = as<std::string>(dynamics_);
  std::string fix = as<std::string>(fix_);
  std::string go_dims = as<std::string>(go_dims_);
  bool immigration = as<bool>(immigration_);
  std::string survey = as<std::string>(survey_);  // andy added this
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
  std::string keyfun = as<std::string>(keyfun_);
  arma::mat a = as<arma::mat>(a_);
  arma::mat u = as<arma::mat>(u_);
  Rcpp::NumericVector db(db_);

  // linear predictors
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::mat omv = arma::ones<arma::colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend")) {
    //    omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
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
  arma::mat sigv = exp(Xsig*beta_sig + Xsig_offset);
  sigv.reshape(T, M);
  arma::mat sig = trans(sigv);
  //* 12/22/2015  *//
 //Immigration
  arma::mat iotav = arma::zeros<arma::colvec>(M*(T-1));
  if(immigration) {
    iotav = exp(Xiota*beta_iota + Xiota_offset);
  }
  iotav.reshape(T-1, M);
  arma::mat iota = arma::trans(iotav);
  //* 12/22/2015 *//
  // format matrices as cubes (shouldn't be done in likelihood)
  arma::icube y(M,J,T);
  arma::icube yna(M,J,T);
  for(int q=0; q<(M*J*T); q++) {
    y(q) = ym(q);
    yna(q) = ynam(q);
  }

  // initialize
  double ll=0.0;
  double w=0.0;
  double f0 = 0.0;

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

  double lower=0.0;
  double upper=0.0;
  double result=0.0;

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
  double part2m3 = 0.0;

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
	// important: fix this!
	f0 = Rf_dnorm4(0.0, 0.0, sig(i,t), false);
	g1_t.zeros();
	// loop over possible value of N at time t
	for(int k=0; k<lk; k++) { // note first index
	  part1 = lgamma(k+1);
	  part2m3 = 0.0;
	  cpsum = 0.0;
	  if((k - yt(i,t)) >= 0) {
	    for(int j=0; j<=J; j++) { // note <=J not <J
	      cp = 0.0;
	      if(j<J) {
               	lower = db[j];
		upper = db[j+1];
		w = upper-lower; // FIXME: Should pass as argument
		result = 0.0;
		if(survey == "line") {
		  result = (R::pnorm(upper, 0.0, sig(i,t), true, false) -
			    R::pnorm(lower, 0.0, sig(i,t), true, false)) / f0;
		  cp = result / w * u(i,j);
		}
		if(survey == "point") {
		  result = sig(i,t)*sig(i,t)*(1-exp(-upper*upper/(2*sig(i,t)*sig(i,t))))-
		    sig(i,t)*sig(i,t)*(1-exp(-lower*lower/(2*sig(i,t)*sig(i,t))));
		  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
		}
	      cp = std::max(cp, DOUBLE_XMIN);
	      cpsum = cpsum+cp;
	      part2 = log(cp) * y(i,j,t);
	      part3 = lgamma(y(i,j,t)+1); // compute outside likelihood
	      part2m3 = part2m3 + (part2 - part3);
	    } else {
	      cp = 1 - cpsum;
	      cp = std::max(cp, DOUBLE_XMIN);
	      part2 = log(cp) * (k - yt(i,t));
	      part3 = lgamma(k-yt(i,t)+1);
	      part2m3 = part2m3 + (part2-part3);
	    }
	  }
	  g1_t(k) = part1 + part2m3;
	  g1_t(k) = exp(g1_t(k));
	}
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
  // important: fix this!
  f0 = Rf_dnorm4(0.0, 0.0, sig(i,first_i), false);

  for(int k=0; k<lk; k++) { // loop over possible values of N
    if((k-yt(i,first_i)) >= 0) {
      part1 = lgamma(k+1); // compute outside likelihood, or ignore
      part2m3 = 0.0;
      cpsum = 0.0;
      for(int j=0; j<=J; j++) {
	cp = 0.0;
	if(j<J) {
          lower = db[j];
	  upper = db[j+1];
	  w = upper-lower; // FIXME: Should pass as argument
	  result = 0.0;
	  if(survey == "line"){
	    result = (R::pnorm(upper, 0.0, sig(i,first_i), true, false) -
		      R::pnorm(lower, 0.0, sig(i,first_i), true, false)) / f0;
	    cp = result / w * u(i,j);
	  }
	  if(survey == "point"){
	    result = sig(i,first_i)*sig(i,first_i) *
	      (1-exp(-upper*upper/(2*sig(i,first_i)*sig(i,first_i)))) -
	      sig(i,first_i)*sig(i,first_i) *
	      (1-exp(-lower*lower/(2*sig(i,first_i)*sig(i,first_i))));
	    cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	  }
	  cp = std::max(cp, DOUBLE_XMIN);
	  //	  if(i==0 && k==(lk-1))
	  //	    Rprintf("        cp : %f \n", cp);
	  cpsum = cpsum+cp;
	  part2 = log(cp) * y(i,j,first_i);
	  part3 = lgamma(y(i,j,first_i)+1); // compute outside likelihood
	  part2m3 = part2m3 + (part2 - part3);
	} else {
	  cp = 1 - cpsum;
	  cp = std::max(cp, DOUBLE_XMIN);
	  //	  if(i==0 && k==(lk-1))
	  //	    Rprintf("        cp : %f \n", cp);
	  part2 = log(cp) * (k - yt(i,first_i));
	  part3 = lgamma(k-yt(i,first_i)+1);
	  part2m3 = part2m3 + (part2 - part3);
	}
      }
      g1(k) = part1 + part2m3;
      g1(k) = exp(g1(k));
    }
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
  //  ll += log(ll_i + DOUBLE_XMIN);
  ll += log(ll_i);
  }

  // Rprintf("        sig(1,1) : %f \n", sig(0,0));
  // Rprintf("        sig(2,2) : %f \n", sig(1,1));
  // Rprintf("        Xsig(1) : %f \n", Xsig(0,0));
  // Rprintf("        Xsig_offset(1) : %f \n", Xsig_offset(0));
  // Rprintf("        beta_sig(1) : %f \n", beta_sig(0));
  // Rprintf("        sigv(1) : %f \n", sigv(0,0));
  // Rprintf("        lower : %f \n", lower);
  // Rprintf("        upper : %f \n", upper);
  // Rprintf("        ier : %d \n", ier);
  // Rprintf("        part1 : %f \n", part1);
  // Rprintf("        part2 : %f \n", part2);
  // Rprintf("        part3 : %f \n", part3);
  // Rprintf("        part2m3 : %f \n", part2m3);
  // Rprintf("        result : %f \n", result);
  // Rprintf("        cpsum : %f \n", cpsum);
  // Rprintf("        cp : %f \n", cp);
  // Rprintf("    Log-likelihood : %f \n", ll);

  return wrap(-ll);
}
