#ifndef _UNMARKED_TRANPROBS_H
#define _UNMARKED_TRANPROBS_H

#include <RcppArmadillo.h>




// constant model
void tp1(arma::mat& g3, int nrI, int nrI1, Rcpp::IntegerVector N, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip, double gam, double om);

// autoregressive + immigration model
void tp2(arma::mat& g3, int lk, double gam, double om, double imm);
 
// trend + immigration model 
void tp3(arma::mat& g3, int lk, double gam, double imm);

// Ricker + immigration model
void tp4(arma::mat& g3, int lk, double gam, double om, double imm);

// Gompertz + immigration model
void tp5(arma::mat& g3, int lk, double gam, double om, double imm);



/*
// constant model
void tp1(arma::mat& g3, int nrI, int nrI1, Rcpp::IntegerVector N, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip, double gam, double om); 

// autoregressive + immigration model
void tp2(arma::mat& g3, int lk, double gam, double om, double imm);
 
// trend + immigration model 
void tp3(arma::mat& g3, int lk, double gam, double imm);

// Ricker + immigration model
void tp4(arma::mat& g3, int lk, double gam, double om, double imm);

// Gompertz + immigration model
void tp5(arma::mat& g3, int lk, double gam, double om, double imm);
*/

#endif
