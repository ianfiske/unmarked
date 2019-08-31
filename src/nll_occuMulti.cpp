#include "nll_occuMulti.h"

using namespace Rcpp;
using namespace arma;

SEXP nll_occuMulti( SEXP fStartR, SEXP fStopR, SEXP dmFr, SEXP dmOccR, 
    SEXP betaR, SEXP dmDetR, SEXP dStartR, SEXP dStopR, SEXP yR, SEXP yStartR, 
    SEXP yStopR, SEXP Iy0r, SEXP zR, SEXP fixed0r){
  
  //Inputs
  IntegerVector fStart(fStartR);
  IntegerVector fStop(fStopR);
  
  //if Matrix is a dependency
  sp_mat dmF = as<sp_mat>(dmFr); //already transposed
  
  //if Matrix not a dependency
  //sp_mat dmF( as<mat>(dmFr) );
  
  int nF = dmF.n_rows;
  List dmOcc(dmOccR);

  colvec beta = as<colvec>(betaR);
  LogicalVector fixed0(fixed0r);
  
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
  int index = 0;
  mat f(N, nF);
  for(int i = 0; i < nF; i++){
    if(fixed0(i)){
      f.col(i) = zeros(N);
    } else {
      mat X = as<mat>(dmOcc[index]);
      f.col(i) = X * beta.subvec(fStart[index], fStop[index]);
      index += 1;
    }
  }

  mat psi = exp( f * dmF ); 
  for(unsigned int i = 0; i < psi.n_rows; i++){
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
