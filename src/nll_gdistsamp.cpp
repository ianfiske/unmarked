#include "nll_gdistsamp.h"

using namespace Rcpp;
using namespace arma;

mat invlogit(const mat& inp ){
  return(1 / (1 + exp(-1 * inp)));
}

vec p_halfnorm(const double& sigma, const std::string& survey, 
               const NumericVector& db, const vec& w, const rowvec& a){
  
  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    double f0 = 2 * R::dnorm(0.0, 0.0, sigma, 0);
    IntegerVector idx1 = Range(1,db.size()-1);
    NumericVector db_sub1 = db[idx1];
    vec p1 = Rcpp::pnorm(db_sub1, 0.0, sigma);
    IntegerVector idx2 = Range(0,db.size()-2);
    NumericVector db_sub2 = db[idx2];
    vec p2 = Rcpp::pnorm(db_sub2, 0.0, sigma);
    vec int_ = 2 * (p1 - p2);
    p = int_ / f0 / w;

  } else if(survey == "point"){
    for (int j; j<J; j++){
      double s2 = pow(sigma,2);
      double p1 = 1 - exp(-pow(db[j+1],2) / (2 * s2));
      double p2 = 1 - exp(-pow(db[j],2) / (2 * s2)); 
      double int_ = s2 * p1 - s2 * p2;
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

SEXP nll_gdistsamp(SEXP beta_, SEXP mixture_, SEXP keyfun_, SEXP survey_,
    SEXP Xlam_, SEXP Xlam_offset_, SEXP A_, SEXP Xphi_, SEXP Xphi_offset_, 
    SEXP Xdet_, SEXP Xdet_offset_, SEXP db_, SEXP a_, SEXP u_, SEXP w_,
    SEXP k_, SEXP lfac_k_, SEXP lfac_kmyt_, SEXP kmyt_, 
    SEXP y_, SEXP naflag_, SEXP fin_, 
    SEXP nP_, SEXP nLP_, SEXP nPP_, SEXP nDP_, SEXP nSP_, SEXP nOP_){

  //Inputs
  vec beta = as<vec>(beta_);
  std::string mixture = as<std::string>(mixture_);
  std::string keyfun = as<std::string>(keyfun_);
  std::string survey = as<std::string>(survey_);

  mat Xlam = as<mat>(Xlam_);
  vec Xlam_offset = as<vec>(Xlam_offset_);
  vec A = as<vec>(A_);
  mat Xphi = as<mat>(Xphi_);
  vec Xphi_offset = as<vec>(Xphi_offset_);
  
  mat Xdet = as<mat>(Xdet_);
  vec Xdet_offset = as<vec>(Xdet_offset_);
  
  NumericVector db(db_);
  vec w = as<vec>(w_);
  mat a = as<mat>(a_);
  mat ut = as<mat>(u_);
  mat u = ut.t();

  IntegerVector k(k_);
  vec lfac_k = as<vec>(lfac_k_);
  cube lfac_kmyt = as<cube>(lfac_kmyt_);
  cube kmyt = as<cube>(kmyt_);

  vec y = as<vec>(y_);
  vec naflag = as<vec>(naflag_);
  mat fin = as<mat>(fin_);

  int nP = as<int>(nP_);
  int nLP = as<int>(nLP_);
  int nPP = as<int>(nPP_);
  int nDP = as<int>(nDP_);
  int nSP = as<int>(nSP_);
  int nOP = as<int>(nOP_);

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nLP - 1) ) + Xlam_offset ) % A;

  int T = Xphi.n_rows / M;
  vec phi = ones(M*T);
  if(T > 1){
    phi = invlogit( Xphi * beta.subvec(nLP, (nLP+nPP-1)) + Xphi_offset);
  }
  
  int R = y.size() / M;
  int J = R / T;
  
  vec detParam(M*T);
  if(keyfun != "uniform"){
    detParam = exp( Xdet * beta.subvec((nLP+nPP),(nLP+nPP+nDP-1)) 
                        + Xdet_offset);
  }
  
  int K = k.size();
  int t_ind = 0;
  int y_ind = 0;
  vec ll(M);
  for (int m=0; m<M; m++){
    vec f(K);
    if(mixture == "P"){
      f = dpois(k, lambda(m));
    } else if(mixture == "NB"){
      f = dnbinom_mu(k, exp(beta(nP-1)), lambda(m));
    }

    mat mn = zeros(K,T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      vec na_sub = naflag.subvec(y_ind, y_stop); 

      if( ! all(na_sub) ){

        vec p1 = lfac_kmyt.subcube(span(m),span(t),span());
        vec p2 = y.subvec(y_ind, y_stop);
        
        //calculate p here
        vec p(J);
        if(keyfun == "uniform"){
          p = ones(J); 
        } else if (keyfun == "halfnorm"){
          //detParam is sigma
          p = p_halfnorm(detParam(t_ind), survey, db, w, a.row(m));
        }
        //other keyfuns

        vec p3 = p % u.col(m) * phi(t_ind);
        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        vec p4 = kmyt.subcube(span(m),span(t),span());

        if( any(na_sub)){
          uvec ids = find(na_sub!=1);
          p2 = p2.elem(ids);
          p3 = p3.elem(ids);
        }

        double p5 = 1 - sum(p3);

        mn.col(t) = lfac_k - p1 + sum(p2 % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      y_ind += J;
    }

    vec g(K);
    for (int i=0; i<K; i++){
      g(i) = exp(sum(mn.row(i)));
      if(!fin(m,i)){
        f(i) = 0;
        g(i) = 0;
      }
    }

    ll(m) = log(sum(f % g));
  }

  return(wrap(-sum(ll)));

}
