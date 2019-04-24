#include "nll_occuMS.h"

using namespace Rcpp;
using namespace arma;

//Build parameter from design matrix and beta
mat get_param(const List& dm_list, const vec& beta, const mat& ind){

  const int m = dm_list.size();
  const mat sz = as<mat>(dm_list[0]); //not ideal
  const int l = sz.n_rows;
  mat out(l,m);

  for(int i=0; i < m; i++){
    mat X = as<mat>(dm_list[i]);
    vec linpred = X * beta.subvec(ind(i,0), ind(i,1));
    out.col(i) = 1 / ( 1 + exp(-linpred));
  }
  
  return(out);
}

//calculate final psi depending on parameterization
mat get_psi(const mat& psi_raw, const std::string& prm){
  int R = psi_raw.n_rows;
  
  if(prm == "multinomial"){
    //[ 1 - psi_1:psi_m, psi_1:psi_m ]
    colvec col1(R);
    for(int i=0; i<R; i++){
      rowvec tmp = psi_raw.row(i);
      col1(i) = 1 - sum(tmp);
    }

    mat out = join_rows(col1, psi_raw);
    return(out);

  } else if(prm == "condbinom"){
    //[ 1-psi, psi * (1-R), psi * R ]
    mat out(R,3);
    for(int i=0; i<R; i++){
      out(i,0) = 1 - psi_raw(i,0);
      out(i,1) = psi_raw(i,0) * (1 - psi_raw(i,1));
      out(i,2) = psi_raw(i,0) * psi_raw(i,1);
    }
    return(out);
  }
}

mat get_sdp(int S, const rowvec& probs, const mat& guide, 
    const std::string& prm){

  mat out = zeros(S,S);

  if(prm == "multinomial"){
    for (int i=0; i<probs.size(); i++){
      out( guide(i,0), guide(i,1) ) = probs(i);
    }
    for (int s=0; s<S; s++){
      out(s,0) = 1 - sum(out.row(s));
    }
    return(out);

  } else if(prm == "condbinom"){
    //input probs order is p_1, p_2, delta
    out(0,0) = 1;
    out(1,0) = 1 - probs(0);
    out(1,1) = probs(0);
    out(2,0) = 1 - probs(1);
    out(2,1) = probs(1) * (1 - probs(2));
    out(2,2) = probs(1) * probs(2);
    return(out);
  }
}

vec get_ph(const int S, const rowvec& y, const mat& probs, 
    const rowvec& navec, const mat& guide, const std::string prm){
  
  int J = probs.n_rows;

  vec out = ones(S);
  for(int j=0; j<J; j++){
    if(! navec(j)){
      rowvec prsub = probs.row(j);
      mat sdp = get_sdp(S, prsub, guide, prm);
      for(int s=0; s<S; s++){
        out(s) = out(s) * sdp(s, y(j));
      }
    }
  }
  return(out);
}

SEXP nll_occuMS( SEXP beta_, SEXP y_, SEXP dm_state_, SEXP dm_det_, 
    SEXP sind_, SEXP dind_, SEXP prm_, SEXP S_, SEXP J_, SEXP N_,
    SEXP naflag_, SEXP guide_){

  //Inputs
  const vec beta = as<vec>(beta_);
  const mat y = as<mat>(y_);
  const List dm_state(dm_state_);
  const List dm_det(dm_det_);

  const mat sind = as<mat>(sind_);
  const mat dind = as<mat>(dind_);
  const mat guide = as<mat>(guide_);

  const std::string prm = as<std::string>(prm_);
  
  const mat naflag = as<mat>(naflag_);
  
  int N = as<int>(N_);
  int S = as<int>(S_);
  int J = as<int>(J_);


  //Get psi values
  const mat raw_psi = get_param(dm_state, beta, sind);
  const mat psi = get_psi(raw_psi, prm);

  //Get p values
  const mat p = get_param(dm_det, beta, dind);

  vec lik = zeros(N);
  int pstart = 0;
  int pend;

  for(int n=0; n<N; n++){
    pend = pstart + J - 1;
    rowvec ysub = y.row(n);
    rowvec nasub = naflag.row(n);
    mat psub = p.rows(span(pstart, pend));
    vec ph = get_ph(S, ysub, psub, nasub, guide, prm);
    lik(n) = dot(psi.row(n), ph);
    pstart += J;
  }

  return(wrap(-sum(log(lik))));

}
