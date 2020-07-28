#include "nll_gdistsamp.h"
#include "distprob.h"

using namespace Rcpp;
using namespace arma;

mat invlogit(const mat& inp ){
  return(1 / (1 + exp(-1 * inp)));
}

SEXP nll_gdistsamp(SEXP beta_, SEXP mixture_, SEXP keyfun_, SEXP survey_,
    SEXP Xlam_, SEXP Xlam_offset_, SEXP A_, SEXP Xphi_, SEXP Xphi_offset_, 
    SEXP Xdet_, SEXP Xdet_offset_, SEXP db_, SEXP a_, SEXP u_, SEXP w_,
    SEXP k_, SEXP lfac_k_, SEXP lfac_kmyt_, SEXP kmyt_, 
    SEXP y_, SEXP naflag_, SEXP fin_, 
    SEXP nP_, SEXP nLP_, SEXP nPP_, SEXP nDP_, SEXP rel_tol_){

  //Inputs
  vec beta = as<vec>(beta_);
  std::string mixture = as<std::string>(mixture_);
  std::string keyfun = as<std::string>(keyfun_);
  std::string survey = as<std::string>(survey_);

  mat Xlam = as<mat>(Xlam_);
  vec Xlam_offset = as<vec>(Xlam_offset_);
  vec A = as<vec>(A_);
  mat Xphi = as<mat>(Xphi_);
  vec Xphi_offset = as<vec>(Xphi_offset_);
  
  mat Xdet = as<mat>(Xdet_);
  vec Xdet_offset = as<vec>(Xdet_offset_);
  
  vec db = as<vec>(db_);
  vec w = as<vec>(w_);
  mat a = as<mat>(a_);
  mat ut = as<mat>(u_);
  mat u = ut.t();

  //IntegerVector k(k_);
  vec k = as<vec>(k_);
  vec lfac_k = as<vec>(lfac_k_);
  cube lfac_kmyt = as<cube>(lfac_kmyt_);
  cube kmyt = as<cube>(kmyt_);

  vec y = as<vec>(y_);
  vec naflag = as<vec>(naflag_);
  mat fin = as<mat>(fin_);

  int nP = as<int>(nP_);
  int nLP = as<int>(nLP_);
  int nPP = as<int>(nPP_);
  int nDP = as<int>(nDP_);
  
  //Integration tol currently unused
  //double rel_tol = as<double>(rel_tol_);

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nLP - 1) ) + Xlam_offset ) % A;

  int T = Xphi.n_rows / M;
  vec phi = ones(M*T);
  if(T > 1){
    phi = invlogit( Xphi * beta.subvec(nLP, (nLP+nPP-1)) + Xphi_offset);
  }
  
  int R = y.size() / M;
  int J = R / T;
  
  vec det_param(M*T);
  if(keyfun != "uniform"){
    det_param = exp( Xdet * beta.subvec((nLP+nPP),(nLP+nPP+nDP-1)) + Xdet_offset);
  }
  
  int K = k.size();
  int t_ind = 0;
  int y_ind = 0;
  vec ll(M);
  for (int m=0; m<M; m++){
    vec f(K);
    if(mixture == "P"){
      for(int l=0; l<K; l++){
        f(l) = R::dpois(k(l), lambda(m), 0);
      }
    } else if(mixture == "NB"){
      double sz = exp(beta(nP-1));
      for(int l=0; l<K; l++){
        f(l) = R::dnbinom_mu(k(l), sz, lambda(m), 0);
      }
    }

    mat mn = zeros(K,T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      vec na_sub = naflag.subvec(y_ind, y_stop); 

      if( ! any(na_sub) ){

        vec p1 = lfac_kmyt.subcube(span(m),span(t),span());
        vec p2 = y.subvec(y_ind, y_stop);
        
        //calculate p
        double scale = 0.0;
        if(keyfun =="hazard"){
          scale = exp(beta(nLP+nPP+nDP));
        }
        vec p = distprob(keyfun, det_param(t_ind), scale, survey, db, 
                          w, a.row(m));
        vec p3 = p % u.col(m) * phi(t_ind);
        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        vec p4 = kmyt.subcube(span(m),span(t),span());

        double p5 = 1 - sum(p3);

        mn.col(t) = lfac_k - p1 + sum(p2 % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      y_ind += J;
    }

    vec g(K);
    for (int i=0; i<K; i++){
      g(i) = exp(sum(mn.row(i)));
      if(!fin(m,i)){
        f(i) = 0;
        g(i) = 0;
      }
    }

    ll(m) = log(sum(f % g) + DOUBLE_XMIN);
  }

  return(wrap(-sum(ll)));

}
