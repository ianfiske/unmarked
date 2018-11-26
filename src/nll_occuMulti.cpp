#include "nll_occuMulti.h"

using namespace Rcpp;
using namespace arma;

SEXP nll_occuMulti( SEXP fStartR, SEXP fStopR, SEXP dmFr, SEXP dmOccR, 
    SEXP betaR, SEXP dmDetR, SEXP dStartR, SEXP dStopR, SEXP yR, SEXP yStartR, 
    SEXP yStopR, SEXP Iy0r, SEXP zR, SEXP fixed0r){
  
  //Inputs
  IntegerVector fStart(fStartR);
  IntegerVector fStop(fStopR);
  mat dmF = as<mat>(dmFr);
  List dmOcc(dmOccR);
  int nF = dmOcc.size();
  colvec beta_tmp = as<colvec>(betaR);
  LogicalVector fixed0(fixed0r);
  int nP = fixed0.size();
  colvec beta(nP);

  int index = 0;
  for(int i = 0; i < nP; i++){
    if(fixed0(i)){
      beta(i) = 0;
    } else {
      beta(i) = beta_tmp(index);
      index += 1;
    }
  }
 
  List dmDet(dmDetR);
  IntegerVector dStart(dStartR);
  IntegerVector dStop(dStopR);
  
  mat y = as<mat>(yR);
  int S = y.n_cols;
  int J = y.n_rows;
  IntegerVector yStart(yStartR);
  IntegerVector yStop(yStopR);
  int N = yStart.size();
  mat Iy0 = as<mat>(Iy0r);

  mat z = as<mat>(zR);
  
  //psi calculation
  mat f(N, nF);
  for(int i = 0; i < nF; i++){
    mat X = as<mat>(dmOcc[i]);
    f.col(i) = X * beta.subvec(fStart[i], fStop[i]);
  }

  mat psi = exp( f * dmF.t() );
  for(int i = 0; i < psi.n_rows; i++){
    psi.row(i) = psi.row(i) / sum( psi.row(i) );
  }

  //p calculation
  mat p(J, S);
  for(int i = 0; i < S; i++){
    mat X = as<mat>(dmDet[i]);
    colvec beta_sub = beta.subvec(dStart[i], dStop[i]);
    p.col(i) = 1.0 / ( 1.0 + exp( -1.0 * (X * beta_sub) ));
  }

  //Probability of detection history
  int M = z.n_rows;
  vec logLik(N);

  for(int i = 0; i < N; i++){
    
    mat ysub = y.rows(yStart[i], yStop[i]);
    mat psub = p.rows(yStart[i], yStop[i]);
    rowvec cdp(S);
    for(int j = 0; j < S; j++){
     cdp(j) = exp( sum( ysub.col(j) % log(psub.col(j)) ) + 
                   sum( (1 - ysub.col(j)) % log( 1 - psub.col(j)) ) );
    }
    
    rowvec prdProbY(M);
    for(int j = 0; j < M; j++){
      prdProbY(j) = prod( z.row(j) % cdp + (1 - z.row(j)) % Iy0.row(i) );
    }

    logLik(i) = log( sum( psi.row(i) % prdProbY ) );

  }

  return(wrap(-1.0 * sum(logLik)));
}
