#include "nll_gdistsamp.h"

using namespace Rcpp;
using namespace arma;

mat invlogit(const mat& inp ){
  return(1 / (1 + exp(-1 * inp)));
}

vec p_halfnorm(const double& sigma, const std::string& survey, 
               const vec& db, const vec& w, const rowvec& a){
  
  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    double f0 = 2 * R::dnorm(0.0, 0.0, sigma, 0);
    int L = db.size();
    vec p1(L-1);
    vec p2(L-1);
    for(int l=1; l<L; l++){
      p1(l-1) = R::pnorm(db(l), 0.0, sigma, 1, 0);
      p2(l-1) = R::pnorm(db(l-1), 0.0, sigma, 1, 0);
    }
    vec int_ = 2 * (p1 - p2);
    p = int_ / f0 / w;

  } else if(survey == "point"){
    for (int j=0; j<J; j++){
      double s2 = pow(sigma,2);
      double p1 = 1 - exp(-pow(db(j+1),2) / (2 * s2));
      double p2 = 1 - exp(-pow(db(j),2) / (2 * s2)); 
      double int_ = s2 * p1 - s2 * p2;
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

vec p_exp(const double& rate, const std::string& survey, const vec& db, 
          const vec& w, const rowvec& a, double& rel_tol){

  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    for(int j=0; j<J; j++){
      double int_ = rate*(1 - exp(-db(j+1)/rate)) - rate*(1-exp(-db(j)/rate));
      p(j) = int_ / w(j);
    }

  } else if(survey == "point"){
    DetExp f(rate, 1);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

vec p_hazard(const double& shape, const double& scale, const std::string& survey, 
          const vec& db, const vec& w, const rowvec& a, double& rel_tol){

  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    DetHaz f(shape, scale, 0);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ / w(j);
    }

  } else if(survey == "point"){
    DetHaz f(shape, scale, 1);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
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
    SEXP nP_, SEXP nLP_, SEXP nPP_, SEXP nDP_, SEXP nSP_, SEXP nOP_,
    SEXP rel_tol_){

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
  
  vec db = as<vec>(db_);
  vec w = as<vec>(w_);
  mat a = as<mat>(a_);
  mat ut = as<mat>(u_);
  mat u = ut.t();

  //IntegerVector k(k_);
  vec k = as<vec>(k_);
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
  
  //Integration tol
  double rel_tol = as<double>(rel_tol_);

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nLP - 1) ) + Xlam_offset ) % A;

  int T = Xphi.n_rows / M;
  vec phi = ones(M*T);
  if(T > 1){
    phi = invlogit( Xphi * beta.subvec(nLP, (nLP+nPP-1)) + Xphi_offset);
  }
  
  int R = y.size() / M;
  int J = R / T;
  
  vec det_param(M*T);
  if(keyfun != "uniform"){
    det_param = exp( Xdet * beta.subvec((nLP+nPP),(nLP+nPP+nDP-1)) + Xdet_offset);
  }
  
  int K = k.size();
  int t_ind = 0;
  int y_ind = 0;
  vec ll(M);
  for (int m=0; m<M; m++){
    vec f(K);
    if(mixture == "P"){
      for(int l=0; l<K; l++){
        f(l) = R::dpois(k(l), lambda(m), 0);
      }
    } else if(mixture == "NB"){
      double sz = exp(beta(nP-1));
      for(int l=0; l<K; l++){
        f(l) = R::dnbinom_mu(k(l), sz, lambda(m), 0);
      }
    }

    mat mn = zeros(K,T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      vec na_sub = naflag.subvec(y_ind, y_stop); 

      if( ! any(na_sub) ){

        vec p1 = lfac_kmyt.subcube(span(m),span(t),span());
        vec p2 = y.subvec(y_ind, y_stop);
        
        //calculate p
        vec p(J);
        if(keyfun == "uniform"){
          p = ones(J); 
        } else if (keyfun == "halfnorm"){
          //det_param is sigma
          p = p_halfnorm(det_param(t_ind), survey, db, w, a.row(m));
        } else if (keyfun == "exp"){
          //det_param is rate
          p = p_exp(det_param(t_ind), survey, db, w, a.row(m), rel_tol);
        } else if (keyfun == "hazard"){
          //det_param is shape
          double scale = exp(beta(nLP+nPP+nDP));
          p = p_hazard(det_param(t_ind), scale, survey, db, w, a.row(m), rel_tol);
        }
        //other keyfuns

        vec p3 = p % u.col(m) * phi(t_ind);
        //the following line causes a segfault only in R CMD check,
        //when kmyt contains NA values
        vec p4 = kmyt.subcube(span(m),span(t),span());

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

    ll(m) = log(sum(f % g) + DOUBLE_XMIN);
  }

  return(wrap(-sum(ll)));

}
