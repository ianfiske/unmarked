#include <RcppArmadillo.h>
#include "distprob.h"
#include "pifun.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_distRemoval(arma::vec beta, arma::uvec n_param, arma::vec yDistance,
    arma::vec yRemoval, int mixture, std::string keyfun,
    arma::mat Xlam, arma::vec A, arma::mat Xphi, arma::mat Xrem,
    arma::mat Xdist, arma::vec db, arma::mat a, arma::mat u, arma::vec w,
    arma::vec k, arma::vec lfac_k,
    arma::cube lfac_kmyt_dist, arma::cube kmyt_dist,
    arma::cube lfac_kmyt_rem, arma::cube kmyt_rem,
    arma::uvec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = Xlam.n_rows;
  int T = Xphi.n_rows / M;
  int Rdist = yDistance.size() / M;
  unsigned Jdist = Rdist / T;
  int Rrem = yRemoval.size() / M;
  unsigned Jrem = Rrem / T;
  int K = k.size() - 1;

  //Abundance
  const vec lambda = exp(Xlam * beta_sub(beta, n_param, 0)) % A;
  double log_alpha = beta_sub(beta, n_param, 1)(0); //length 1 vector

  //Availability
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(Xphi * beta_sub(beta, n_param, 2));
  }

  //Distance sampling detection
  vec dist_param(M*T);
  if(keyfun != "uniform"){
    dist_param = exp(Xdist * beta_sub(beta, n_param, 3));
  }
  double scale = exp(beta_sub(beta, n_param, 4)(0));

  //Removal sampling detection
  vec remP(M * T * Jrem);
  remP = inv_logit(Xrem * beta_sub(beta, n_param, 5));

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){

    int t_ind = i * T;
    int yd_ind = i * T * Jdist;
    int yr_ind = i * T * Jrem;

    vec yd_sub(Jdist);
    vec yr_sub(Jrem);

    mat mn_dist = zeros(K+1, T);
    mat mn_rem = zeros(K+1, T);

    for(int t=0; t<T; t++){
      int yd_stop = yd_ind + Jdist - 1;
      yd_sub = yDistance.subvec(yd_ind, yd_stop);

      int yr_stop = yr_ind + Jrem - 1;
      yr_sub = yRemoval.subvec(yr_ind, yr_stop);

      // Distance sampling
      uvec not_missing = find_finite(yd_sub);

      vec p1(K+1);
      vec p4(K+1);

      if(not_missing.size() == Jdist){
        p1 = lfac_kmyt_dist.subcube(span(i),span(t),span());

        vec p = distprob(keyfun, dist_param(t_ind), scale, "point", db,
                          w, a.row(i));
        vec p3 = p % u.col(i) * phi(t_ind);

        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        p4 = kmyt_dist.subcube(span(i),span(t),span());

        double p5 = 1 - sum(p3);
        mn_dist.col(t) = lfac_k - p1 + sum(yd_sub % log(p3)) + p4 * log(p5);
      }

      // Removal sampling
      not_missing = find_finite(yr_sub);
      if(not_missing.size() == Jrem){
        p1 = lfac_kmyt_rem.subcube(span(i),span(t),span());
        vec p3 = piFun( remP.subvec(yr_ind, yr_stop), "removalPiFun") * phi(t_ind);
        //segfault on this line
        //vec p4 = kmyt_rem.subcube(span(i),span(t),span());
        p4 = kmyt_rem.slice(i).col(t);
        double p5 = 1 - sum(p3);

        mn_rem.col(t) = lfac_k - p1 + sum(yr_sub % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      yd_ind += Jdist;
      yr_ind += Jrem;
    }

    double site_lp = 0.0;
    for (int j=Kmin(i); j<(K+1); j++){
      site_lp += N_density(mixture, j, lambda(i), log_alpha) *
        exp(sum(mn_dist.row(j))) *
        exp(sum(mn_rem.row(j)));
    }

    loglik += log(site_lp + DOUBLE_XMIN);

  }

  return -loglik;

}
