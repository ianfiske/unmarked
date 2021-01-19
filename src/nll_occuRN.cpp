#include <RcppArmadillo.h>
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double lp_site_occuRN(const rowvec y, double lam, const vec q, int K, int Kmin){

  //y must be arma::vec for this to work (not uvec/ivec)
  const uvec fin = find_finite(y);
  if(fin.size() == 0) return 0.0;

  double f, g, p, out = 0.0;

  for (int k=Kmin; k<(K+1); k++){
    f = Rf_dpois(k, lam, 0);
    g = 0.0;
    for (unsigned j=0; j<fin.size(); j++){
      p = 1 - pow(q(fin(j)), k);
      g += Rf_dbinom(y(fin(j)), 1, p, true);
    }
    out += f * exp(g);
  }
  return log(out + DOUBLE_XMIN);
}

// [[Rcpp::export]]
double nll_occuRN(const arma::vec beta, const arma::uvec n_param, const arma::mat y,
                  const arma::mat X, const arma::mat V, const arma::vec X_offset,
                  const arma::vec V_offset, int K, const arma::uvec Kmin,
                  int threads){

  int M = y.n_rows;
  int J = y.n_cols;

  const vec lam = exp(X*beta_sub(beta, n_param, 0) + X_offset);
  const vec q = 1 - inv_logit(V*beta_sub(beta, n_param, 1) + V_offset);

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){
    int pstart = i * J;
    int pstop = i * J + J - 1;
    loglik += lp_site_occuRN(y.row(i), lam(i), q.subvec(pstart, pstop),
                             K, Kmin(i));
  }

  return -loglik;

}
