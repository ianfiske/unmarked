#include <RcppArmadillo.h>
#include "distr.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp ;
using namespace arma ;

// [[Rcpp::export]]
double nll_gpcount(arma::mat ym, arma::mat Xlam, arma::mat Xphi, arma::mat Xp,
    arma::vec beta_lam, arma::vec beta_phi, arma::vec beta_p, double log_alpha,
    arma::vec Xlam_offset, arma::vec Xphi_offset, arma::vec Xp_offset,
    int M, std::string mixture, int T, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  arma::mat yimax = max(ym, 1);
  int lM = M+1;
  double alpha = 0.0;
  if(mixture=="NB")
    alpha = exp(log_alpha);
  else if(mixture=="ZIP")
    alpha = 1.0/(1.0+exp(-log_alpha));
  int R = ym.n_rows;
  int J = ym.n_cols / T;
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec logit_phi = Xphi*beta_phi + Xphi_offset;
  arma::mat phiv = 1.0/(1.0+exp(-logit_phi));
  arma::colvec logit_p = Xp*beta_p + Xp_offset;
  arma::mat pv = 1.0/(1.0+exp(-logit_p));
  phiv.reshape(T,R);
  arma::mat phi = trans(phiv);
  pv.reshape(J*T,R);
  arma::mat pm = trans(pv);
  arma::cube y(R,J,T);
  arma::cube p(R,J,T);
  for(int q=0; q<(R*J*T); q++) {
    y(q) = ym(q);
    p(q) = pm(q);
  }
  //  Rprintf("made it %f \\n", 1.);
  double loglik=0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for(int i=0; i<R; i++) {

    double g=0.0, h=0.0;
    arma::vec f = arma::zeros<arma::vec>(lM);
    arma::mat gh = arma::zeros<arma::mat>(lM, lM);
    arma::vec ghi = arma::zeros<arma::vec>(lM);

    //    Rprintf("log-like %f \\n", L);
    for(int m=0; m<lM; m++) {
      if(m < yimax(i))
	f(m) = log(0.0);
      else if(mixture=="P")
	f(m) = Rf_dpois(m, lam(i), true);
      else if(mixture=="NB")
        f(m) = dnbinom_mu(m, alpha, lam(i), true);
      else if(mixture=="ZIP")
	f(m) = log(dzip(m, lam(i), alpha));
    }
    ghi.zeros();
    for(int t=0; t<T; t++) {
      gh.zeros();
      for(int m=0; m<lM; m++) { // FIXME: increment from max(y[i,])
	for(int n=0; n<lM; n++) {
	  if((n > m) || (m < yimax(i))) {
	    gh(n,m) = log(0.0);
	    continue;
	  }
	  g = 0.0;
	  if(arma::is_finite(phi(i,t)))
	    g = Rf_dbinom(n, m, phi(i,t), true);
	  h = 0.0;
	  for(int j=0; j<J; j++) {
	    if(arma::is_finite(y(i,j,t))) { // true if not NA, NaN, +/-Inf
	      h += Rf_dbinom(y(i,j,t), n, p(i,j,t), true);
	    }
	  }
	  gh(n,m) = g + h;
	}
	ghi(m) += log(arma::accu(exp(gh.col(m)))); // sum over N(t)
      }
    }
    loglik += log(arma::accu(exp(f + ghi))); // sum over M
  }
  return -loglik;
}
