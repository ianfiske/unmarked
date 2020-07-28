#include "nll_gmultmix.h"
#include "pifun.h"

using namespace Rcpp;
using namespace arma;

//mat inv_logit( mat inp ){
//  return(1 / (1 + exp(-1 * inp)));
//}


SEXP nll_gmultmix(SEXP betaR, SEXP mixtureR, SEXP pi_funR, 
    SEXP XlamR, SEXP Xlam_offsetR, SEXP XphiR, SEXP Xphi_offsetR, SEXP XdetR, 
    SEXP Xdet_offsetR, SEXP kR, SEXP lfac_kR, SEXP lfac_kmytR, SEXP kmytR, 
    SEXP yR, SEXP naflagR, SEXP finR, SEXP nPr, SEXP nLPr, SEXP nPPr, SEXP nDPr){

  //Inputs
  vec beta = as<vec>(betaR);
  std::string mixture = as<std::string>(mixtureR);
  std::string pi_fun = as<std::string>(pi_funR);

  mat Xlam = as<mat>(XlamR);
  vec Xlam_offset = as<vec>(Xlam_offsetR);
  mat Xphi = as<mat>(XphiR);
  vec Xphi_offset = as<vec>(Xphi_offsetR);
  mat Xdet = as<mat>(XdetR);
  vec Xdet_offset = as<vec>(Xdet_offsetR);

  IntegerVector k(kR);
  vec lfac_k = as<vec>(lfac_kR);
  cube lfac_kmyt = as<cube>(lfac_kmytR);
  cube kmyt = as<cube>(kmytR);

  vec y = as<vec>(yR);
  vec naflag = as<vec>(naflagR);
  mat fin = as<mat>(finR);

  int nP = as<int>(nPr);
  int nLP = as<int>(nLPr);
  int nPP = as<int>(nPPr);
  int nDP = as<int>(nDPr);

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
  int t_ind = 0;
  int y_ind = 0;
  int p_ind = 0;
  vec ll(M);
  for (int m=0; m<M; m++){
    vec f(K);
    if(mixture == "P"){
      f = dpois(k, lambda(m));
    } else if(mixture == "NB"){
      f = dnbinom_mu(k, exp(beta(nP-1)), lambda(m));
    }

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

    ll(m) = log(sum(f % g));
  }

  return(wrap(-sum(ll)));

}
