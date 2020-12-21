#include <RcppArmadillo.h>
#include "pifun.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gmultmix(arma::vec beta, arma::uvec n_param, arma::vec y,
                    int mixture, std::string pi_fun, arma::mat Xlam,
                    arma::vec Xlam_offset, arma::mat Xphi, arma::vec Xphi_offset,
                    arma::mat Xdet, arma::vec Xdet_offset,
                    arma::vec k, arma::vec lfac_k, arma::cube lfac_kmyt,
                    arma::cube kmyt, arma::vec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = Xlam.n_rows;
  int T = Xphi.n_rows / M;
  int J = Xdet.n_rows / (M * T);
  int R = y.size() / (M * T);
  int K = k.size() - 1;

  vec lambda = exp(Xlam * beta_sub(beta, n_param, 0) + Xlam_offset);
  double log_alpha = beta_sub(beta, n_param, 3)(0);

  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(Xphi * beta_sub(beta, n_param, 1) + Xphi_offset);
  }

  vec p = inv_logit(Xdet * beta_sub(beta, n_param, 2) + Xdet_offset);

  double loglik = 0.0;
  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){

    int y_ind = i * T * R;
    int p_ind = i * T * J;
    int t_ind = i * T;

    mat A = zeros(K+1,T);
    vec y_sub;
    uvec fin;
    for(int t=0; t<T; t++){
      int y_stop = y_ind + R - 1;
      int p_stop = p_ind + J - 1;

      y_sub = y.subvec(y_ind, y_stop);
      fin = find_finite(y_sub);
      if(fin.size() == 0) continue;

      vec p1 = lfac_kmyt.subcube(span(i),span(t),span());
      vec p3 = piFun( p.subvec(p_ind, p_stop), pi_fun ) * phi(t_ind);
      //the following line causes a segfault only in R CMD check,
      //when kmyt contains NA values
      vec p4 = kmyt.subcube(span(i),span(t),span());

      y_sub = y_sub.elem(fin);
      p3 = p3.elem(fin);
      double p5 = 1 - sum(p3);

      A.col(t) = lfac_k - p1 + sum(y_sub % log(p3)) + p4 * log(p5);

      t_ind += 1;
      y_ind += R;
      p_ind += J;
    }

    double site_lp = 0.0;
    for (int j=Kmin(i); j<(K+1); j++){
      site_lp += N_density(mixture, j, lambda(i), log_alpha) *
        exp(sum(A.row(j)));
    }

    loglik += log(site_lp);
  }

  return(-loglik);

}
