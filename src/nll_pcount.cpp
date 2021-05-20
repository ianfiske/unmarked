#include <RcppArmadillo.h>
#include <float.h>
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double lp_site_pcount(const rowvec y, int mixture, double lam, double log_alpha,
                      const vec p, int K, int Kmin){

  //y must be arma::vec for this to work (not uvec/ivec)
  uvec fin = find_finite(y);
  if(fin.size() == 0) return 0.0;

  double f, g, out = 0.0;
  for (int k=Kmin; k<(K+1); k++){
    f = N_density(mixture, k, lam, log_alpha);
    g = 0.0;
    for (unsigned j=0; j<fin.size(); j++){
      g += Rf_dbinom(y(fin(j)), k, p(fin(j)), true);
    }
    out += f * exp(g);
  }
  return log(out + DBL_MIN);
}


// [[Rcpp::export]]
double nll_pcount(const arma::vec beta, const arma::uvec n_param, const arma::mat y,
                  const arma::mat X, const arma::mat V, const arma::vec X_offset,
                  const arma::vec V_offset, int K, const arma::uvec Kmin,
                  int mixture, int threads){

  int M = y.n_rows;
  int J = y.n_cols;

  vec lam = exp(X*beta_sub(beta, n_param, 0) + X_offset);
  vec p = inv_logit(V*beta_sub(beta, n_param, 1) + V_offset);
  double log_alpha = beta_sub(beta, n_param, 2)(0);

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  double loglik = 0.0;

  //This will compile but throw an unknown pragma warning if no openMP
  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){
    int pstart = i * J;
    int pstop = i * J + J - 1;
    loglik += lp_site_pcount(y.row(i), mixture, lam(i), log_alpha,
                             p.subvec(pstart, pstop), K, Kmin(i));
  }

  return -loglik;

}
