#include "get_mlogit.h"

using namespace Rcpp;
using namespace arma;

SEXP get_mlogit(SEXP lp_mat_, SEXP type_, SEXP S_, SEXP guide_){
  
  const mat lp_mat = as<mat>(lp_mat_);
  const std::string type = as<std::string>(type_);
  int S = as<int>(S_);
  const mat guide = as<mat>(guide_);

  int R = lp_mat.n_rows;
  int C = lp_mat.n_cols;

  if(type == "state"){

    mat out = exp(lp_mat);
    for(int r=0; r<R; r++){
      double row_sum = sum(out.row(r)) + 1;
      for(int c=0; c<C; c++){
        out(r,c) = out(r,c) / row_sum;
      }
    }
    return(wrap(out));

  } else if(type == "det"){
    
    mat out(R,C);
    for(int r=0; r<R; r++){
      mat sdp = zeros(S,S);
      sdp.col(0) = ones(S);
      for(int c=0; c<C; c++){
        sdp( guide(c,0), guide(c,1)) = exp(lp_mat(r,c));
      }
      for(int s=0; s<S; s++){
        double row_sum = sum(sdp.row(s));
        for (int j=0; j<S; j++){
          sdp(s,j) = sdp(s,j) / row_sum;
        }
      } 
      for(int c=0; c<C; c++){
        out(r,c) = sdp( guide(c,0), guide(c,1) );
      }
    }
    return(wrap(out));
  }
}
