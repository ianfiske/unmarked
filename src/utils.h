#ifndef _unmarked_UTILS_H
#define _unmarked_UTILS_H

#include <RcppArmadillo.h>

arma::mat inv_logit( arma::mat inp );

//arma::vec inv_logit( arma::vec inp );

double inv_logit(double x);

arma::vec beta_sub(arma::vec beta, arma::uvec n_param, unsigned idx);

#endif
