#include <RcppArmadillo.h>
#include <omp.h>
#include "pifun.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gmultmix(arma::vec beta, std::string mixture, std::string pi_fun,
                    arma::mat Xlam, arma::vec Xlam_offset, arma::mat Xphi,
                    arma::vec Xphi_offset, arma::mat Xdet, arma::vec Xdet_offset,
                    Rcpp::IntegerVector k, arma::vec lfac_k, arma::cube lfac_kmyt,
                    arma::cube kmyt, arma::vec y, arma::vec naflag, arma::mat fin,
                    int nP, int nLP, int nPP, int nDP, int threads){

  omp_set_num_threads(threads);

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nLP - 1) ) + Xlam_offset );

  int T = Xphi.n_rows / M;
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit( Xphi * beta.subvec(nLP, (nLP+nPP-1)) + Xphi_offset);
  }

  int J = Xdet.n_rows / (M * T);
  int R = y.size() / (M * T);
  vec p = inv_logit( Xdet * beta.subvec((nLP+nPP),(nLP+nPP+nDP-1)) + Xdet_offset);

  int K = k.size();

  //vec ll(M);
  double loglik = 0.0;
  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int m=0; m<M; m++){
    vec f(K);
    if(mixture == "P"){
      f = dpois(k, lambda(m));
    } else if(mixture == "NB"){
      f = dnbinom_mu(k, exp(beta(nP-1)), lambda(m));
    }

    int y_ind = m * T * R;
    int p_ind = m * T * J;
    int t_ind = m * T;

    mat A = zeros(K,T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + R - 1;
      int p_stop = p_ind + J - 1;
      vec na_sub = naflag.subvec(y_ind, y_stop);

      if( ! all(na_sub) ){

        vec p1 = lfac_kmyt.subcube(span(m),span(t),span());
        vec p2 = y.subvec(y_ind, y_stop);
        vec p3 = piFun( p.subvec(p_ind, p_stop), pi_fun ) * phi(t_ind);
        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        vec p4 = kmyt.subcube(span(m),span(t),span());

        if( any(na_sub)){
          uvec ids = find(na_sub!=1);
          p2 = p2.elem(ids);
          p3 = p3.elem(ids);
        }

        double p5 = 1 - sum(p3);

        A.col(t) = lfac_k - p1 + sum(p2 % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      y_ind += R;
      p_ind += J;
    }

    vec g(K);
    for (int i=0; i<K; i++){
      g(i) = exp(sum(A.row(i)));
      if(!fin(m,i)){
        f(i) = 0;
        g(i) = 0;
      }
    }

    loglik += log(sum(f % g));
  }

  return(-loglik);

}
