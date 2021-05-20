#include "utils.h"

arma::mat inv_logit( arma::mat inp ){
  return(1 / (1 + exp(-1 * inp)));
}

//arma::vec inv_logit( arma::vec inp ){
//  return(1 / (1 + exp(-1 * inp)));
//}

double inv_logit(double x){
  return 1 / (1 + exp(-1 * x));
}

arma::vec beta_sub(arma::vec beta, arma::uvec n_param, unsigned idx){
  unsigned np = n_param(idx);
  if(np == 0) return arma::zeros(1);
  unsigned end = sum(n_param.subvec(0, idx)) - 1;
  unsigned start = end - np + 1;
  return beta.subvec(start, end);
}
