#include <RcppArmadillo.h>
#include "distprob.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gdistsamp(arma::vec beta, arma::uvec n_param, arma::vec y,
    int mixture, std::string keyfun, std::string survey,
    arma::mat Xlam, arma::vec Xlam_offset, arma::vec A, arma::mat Xphi,
    arma::vec Xphi_offset, arma::mat Xdet, arma::vec Xdet_offset, arma::vec db,
    arma::mat a, arma::mat u, arma::vec w, arma::vec k, arma::vec lfac_k,
    arma::cube lfac_kmyt, arma::cube kmyt, arma::uvec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = Xlam.n_rows;
  int T = Xphi.n_rows / M;
  int R = y.size() / M;
  unsigned J = R / T;
  int K = k.size() - 1;

  //Abundance
  const vec lambda = exp(Xlam * beta_sub(beta, n_param, 0) + Xlam_offset) % A;
  double log_alpha = beta_sub(beta, n_param, 4)(0); //length 1 vector

  //Availability
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(Xphi * beta_sub(beta, n_param, 1) + Xphi_offset);
  }

  //Detection
  vec det_param(M*T);
  if(keyfun != "uniform"){
    det_param = exp(Xdet * beta_sub(beta, n_param, 2) + Xdet_offset);
  }
  double scale = exp(beta_sub(beta, n_param, 3)(0));

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){

    int t_ind = i * T;
    int y_ind = i * T * J;

    vec y_sub(J);

    mat mn = zeros(K+1, T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      y_sub = y.subvec(y_ind, y_stop);
      uvec not_missing = find_finite(y_sub);

      if(not_missing.size() == J){

        vec p1 = lfac_kmyt.subcube(span(i),span(t),span());
        vec p = distprob(keyfun, det_param(t_ind), scale, survey, db,
                          w, a.row(i));
        vec p3 = p % u.col(i) * phi(t_ind);
        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        vec p4 = kmyt.subcube(span(i),span(t),span());

        double p5 = 1 - sum(p3);

        mn.col(t) = lfac_k - p1 + sum(y_sub % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      y_ind += J;
    }

    double site_lp = 0.0;
    for (int j=Kmin(i); j<(K+1); j++){
      site_lp += N_density(mixture, j, lambda(i), log_alpha) *
        exp(sum(mn.row(j)));
    }

    loglik += log(site_lp + DOUBLE_XMIN);

  }

  return -loglik;

}
