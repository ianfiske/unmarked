#include "nll_nmixTTD.h"

using namespace Rcpp;
using namespace arma;

SEXP nll_nmixTTD( SEXP beta_, SEXP y_, SEXP delta_,
    SEXP W_, SEXP V_, SEXP pind_, SEXP dind_,
    SEXP tdist_, SEXP N_, SEXP J_, SEXP K_, SEXP naflag_){

  //Inputs
  const vec beta = as<vec>(beta_);
  const vec y = as<vec>(y_);
  const vec delta = as<vec>(delta_);
  const vec naflag = as<vec>(naflag_);
  const mat W = as<mat>(W_);
  const mat V = as<mat>(V_);
  const vec pind = as<vec>(pind_);
  const vec dind = as<vec>(dind_);
  const std::string tdist = as<std::string>(tdist_);

  int N = as<int>(N_);
  int J = as<int>(J_);
  int K = as<int>(K_);

  //Get abundance lambda values
  const vec lamN = exp(W * beta.subvec(pind(0), pind(1)));

  //Get detection lambda values
  const vec lamP = exp(V * beta.subvec(dind(0), dind(1)));

  //Get shape if weibull
  double shp = 1.0;
  if(tdist == "weibull"){
    shp = exp(beta(beta.size() - 1));
  }

  vec lik(N);
  int ystart = 0;
  int yend;
  for (int n=0; n<N; n++){

    vec pK(K+1);
    for (int k=0; k<(K+1); k++){
      pK(k) = R::dpois(k, lamN(n), 0);
    }

    yend = ystart + J - 1;
    vec nasub = naflag.subvec(ystart, yend);
    vec dsub = delta.subvec(ystart, yend);
    vec lamsub = lamP.subvec(ystart, yend);
    vec ysub = y.subvec(ystart, yend);

    if(any(nasub)){
      uvec ids = find(nasub != 1);
      dsub = dsub.elem(ids);
      lamsub = lamsub.elem(ids);
      ysub = ysub.elem(ids);
    }
    int ys = ysub.size();

    vec pY(K+1);
    pY(0) = 1 - any(dsub);

    if(tdist == "weibull"){
      for (int k=1; k<(K+1); k++){
        vec e_lamt(lamsub.size());
        for (int i=0; i<ys; i++){
          double lam = lamsub(i) * k;
          e_lamt(i) = pow(shp*lam*pow(lam*ysub(i),(shp-1)),dsub(i)) *
            exp(-1*pow(lam*ysub(i),shp));
        }
        pY(k) = prod(e_lamt);
      }
    } else {
      for (int k=1; k<(K+1); k++){
        vec e_lamt(lamsub.size());
        for(int i=0; i<ys; i++){
          //exponential
          double lam = lamsub(i) * k;
          e_lamt(i) = pow(lam,dsub(i)) * exp(-lam*ysub(i));
        }
        pY(k) = prod(e_lamt);
      }
    }

    ystart += J;

    lik(n) = dot(pK, pY);

  }

  return(wrap(-sum(log(lik))));

}
