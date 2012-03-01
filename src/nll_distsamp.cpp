#include "nll_distsamp.h"

using namespace Rcpp ;

SEXP nll_distsamp( SEXP y_, SEXP X_, SEXP V_, SEXP beta_lam_, SEXP beta_sig_, SEXP X_offset_, SEXP V_offset_ ) { // also need A,u,a,output

  arma::imat y = as<arma::imat>(y_);
  arma::mat X = as<arma::mat>(X_);
  arma::mat V = as<arma::mat>(V_);
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_sig = as<arma::colvec>(beta_sig_);
  arma::colvec X_offset = as<arma::colvec>(X_offset_);
  arma::colvec V_offset = as<arma::colvec>(V_offset_);
  int R = y.n_rows;
  int J = y.n_cols;
  arma::colvec lam = exp(X*beta_lam + X_offset);
  arma::colvec sigv = exp(V*beta_sig + V_offset);
  // need to add output option to model density instead of abundance
  if(output=="density")
    lam = lam % A; // double check %
  sigv.reshape(J,R);
  arma::mat sig = trans(sigv);
  for(int i=0; i<R; i++) {
    for(int j=0; j<J; j++) {
    }
  }
  return wrap(-ll);
}
