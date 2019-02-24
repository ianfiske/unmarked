#include "nll_occuRN.h"

using namespace Rcpp;
using namespace arma;

mat inv_logitRN( mat inp ){
  return(1 / (1 + exp(-1 * inp)));
}


SEXP nll_occuRN(SEXP betaR, 
    SEXP XlamR, SEXP Xlam_offsetR, SEXP XdetR, SEXP Xdet_offsetR,  
    SEXP KR, SEXP yR, SEXP navecR, SEXP nPr, SEXP nOPr){

  //Inputs
  vec beta = as<vec>(betaR);

  mat Xlam = as<mat>(XlamR);
  vec Xlam_offset = as<vec>(Xlam_offsetR);
  mat Xdet = as<mat>(XdetR);
  vec Xdet_offset = as<vec>(Xdet_offsetR);

  int K = as<int>(KR);
  
  vec y = as<vec>(yR);
  vec navec = as<vec>(navecR);

  int nP = as<int>(nPr);
  int nOP = as<int>(nOPr);

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nOP - 1) ) + Xlam_offset );
  
  int R = y.size();
  int J = R / M;
  vec p = inv_logitRN( Xdet * beta.subvec(nOP,(nP-1)) + Xdet_offset);
  

  int p_ind = 0;
  int p_stop = 0;

  mat p_mat(p.size(),(K+1));
  mat cp_mat = ones(p.size(),(K+1));
  mat lam_in(M,(K+1));
  mat cp_in(M,(K+1));
  vec subcp(J);
  for (int i=0; i<(K+1); i++){

    p_mat.col(i) = 1 - pow((1 - p),i);

    for (int r=0; r<R; r++){
      if(!navec(r)){
        cp_mat(r,i) = pow(p_mat(r,i), y(r)) * pow((1 - p_mat(r,i)),(1-y(r)));
      }
    }

    p_ind = 0;
    for (int m=0; m<M; m++){
      p_stop = p_ind + J - 1;
      lam_in(m,i) = R::dpois(i, lambda(m), 0);
      subcp = cp_mat.submat(span(p_ind,p_stop),span(i));
      cp_in(m,i) = prod(subcp);
      p_ind += J;
    }

  }
  
  mat cp_lam = cp_in % lam_in;
  vec ll(M);
  for (int m=0; m<M; m++){
    ll(m) = log(accu(cp_lam.row(m)));
  }

  return(wrap(-accu(ll)));

}
