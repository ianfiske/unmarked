#include "nll_distsampOpen.h"

using namespace Rcpp;
using namespace arma;

SEXP nll_distsampOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_, 
    SEXP Xsig_, SEXP Xiota_, 
    //SEXP beta_lam_, SEXP beta_gam_, SEXP beta_om_, 
    //SEXP beta_sig_, SEXP beta_iota_, SEXP log_alpha_,
    SEXP beta_, SEXP beta_ind_,
    SEXP Xlam_offset_, 
    SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xsig_offset_, SEXP Xiota_offset_, 
    SEXP ytna_, SEXP lk_, SEXP mixture_, SEXP first_, SEXP last_, 
    SEXP M_, SEXP T_, SEXP delta_, SEXP dynamics_, SEXP survey_, 
    SEXP fix_, SEXP go_dims_, SEXP immigration_, SEXP I_, SEXP I1_, SEXP Ib_, 
    SEXP Ip_, 
    //SEXP scale_, 
    SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_,
    SEXP lfac_k_, SEXP lfac_kmyt_, SEXP kmyt_, SEXP fin_ ) {
  
  //Indices
  int lk = as<int>(lk_);
  Rcpp::IntegerVector N = seq_len(lk)-1;
  int M = as<int>(M_);
  int T = as<int>(T_);
  cube y = as<cube>(y_);
  imat yt = as<imat>(yt_);
  Rcpp::IntegerVector first(first_);
  Rcpp::IntegerVector last(last_);
  arma::imat ytna = as<arma::imat>(ytna_); // y[i,,t] are all NA
  arma::imat delta = as<arma::imat>(delta_);
  vec lfac_k = as<vec>(lfac_k_);
  cube lfac_kmyt = as<cube>(lfac_kmyt_);
  cube kmyt = as<cube>(kmyt_);
  cube fin = as<cube>(fin_);
 
  //Covariate matrices
  mat Xlam = as<mat>(Xlam_);
  mat Xgam = as<mat>(Xgam_);
  mat Xom = as<mat>(Xom_);
  mat Xsig = as<mat>(Xsig_);
  mat Xiota = as<mat>(Xiota_);

  //Model types
  std::string keyfun = as<std::string>(keyfun_);
  std::string mixture = as<std::string>(mixture_);
  std::string dynamics = as<std::string>(dynamics_);
  std::string fix = as<std::string>(fix_);
  std::string go_dims = as<std::string>(go_dims_);
  bool immigration = as<bool>(immigration_);
  std::string survey = as<std::string>(survey_);

  //Parameter vectors
  vec beta = as<vec>(beta_);
  umat bi = as<umat>(beta_ind_);
  vec beta_lam = beta.subvec(bi(0,0), bi(0,1));
  vec beta_gam = beta.subvec(bi(1,0), bi(1,1));
  vec beta_om = beta.subvec(bi(2,0), bi(2,1));
  
  vec beta_sig;
  if(keyfun != "uniform"){
    beta_sig = beta.subvec(bi(3,0), bi(3,1));
  }

  vec beta_iota;
  if(immigration){
    beta_iota = beta.subvec(bi(4,0), bi(4,1));
  }

  double scale;
  if(keyfun == "hazard"){
    scale = exp(beta(bi(5,0)));
  }

  //Get 2nd abundance dist parameter set up if necessary
  //double log_alpha = as<double>(log_alpha_);
  //double alpha=0.0, psi=0.0;
  double alpha, psi;
  if(mixture=="NB"){
    //alpha = exp(log_alpha);
    alpha = exp(beta(bi(6,0)));
  } else if(mixture=="ZIP"){
    //psi = 1.0/(1.0+exp(-log_alpha));
    psi = 1.0 / (1.0 + exp(-1 * beta(bi(6,0))));
  }
 
  colvec Xlam_offset = as<colvec>(Xlam_offset_);
  colvec Xgam_offset = as<colvec>(Xgam_offset_);
  colvec Xom_offset = as<colvec>(Xom_offset_);
  colvec Xsig_offset = as<colvec>(Xsig_offset_);
  colvec Xiota_offset = as<colvec>(Xiota_offset_);  
  

  //Other settings

  imat I = as<arma::imat>(I_);
  imat I1 = as<arma::imat>(I1_);
  List Ib(Ib_);
  List Ip(Ip_);
  int nrI = I.n_rows;
  int nrI1 = I1.n_rows;

  //Distance sampling info
  mat a = as<mat>(a_);
  mat u = as<mat>(u_);
  vec w = as<vec>(w_); 
  vec db = as<vec>(db_);

  //Linear predictors
  //Lambda
  colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  
  //Omega
  mat omv = ones<colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend")) {
    if((dynamics == "ricker")  || (dynamics == "gompertz")){
        omv = exp(Xom*beta_om + Xom_offset);
    } else if((dynamics == "constant")  || (dynamics == "autoreg") || 
              (dynamics == "notrend")){
        omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
    }
  }
  omv.reshape(T-1, M);
  mat om = trans(omv);
  
  //Gamma
  mat gam = zeros(M, T-1);
  if(dynamics == "notrend") {
    mat lamMat = repmat(lam, 1, T-1);
    gam = (1-om) % lamMat;
  } else if (fix != "gamma") {
    mat gamv = exp(Xgam*beta_gam + Xgam_offset);
    gamv.reshape(T-1, M);
    gam = trans(gamv);
  }

  //Sigma (detection param)
  mat sig(M, T);
  if(keyfun != "uniform"){
    mat sigv = exp(Xsig*beta_sig + Xsig_offset);
    sigv.reshape(T, M);
    sig = trans(sigv);
  }

 //Immigration
  mat iotav = zeros<colvec>(M*(T-1));
  if(immigration) {
    iotav = exp(Xiota*beta_iota + Xiota_offset);
  }
  iotav.reshape(T-1, M);
  mat iota = trans(iotav);

  // initialize
  double ll=0.0;

  double ll_i=0.0;
  int first_i=0;
  int last_i=0;
  int first1=0;
  colvec g_star = ones(lk);
  mat g3 = zeros(lk, lk);
  mat g3_d = zeros(lk, lk);
  colvec g1_t = zeros(lk);
  colvec g1_t_star = zeros(lk);
  colvec g1 = zeros(lk);
  colvec g1_star = zeros(lk);
  cube g3_t = zeros(lk,lk,T-1);

  // shouldn't be done in likelihood
  for(int i=0; i<M; i++) {
    if(first[i]==1) {
      first1=i; // site with no missing values at t=1
      break;
    }
  }

  // compute g3 if there are no covariates of omega/gamma
  if(go_dims == "scalar") {
    if(dynamics=="constant" || dynamics=="notrend")
      tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,0), om(first1,0));
    else if(dynamics=="autoreg")
      tp2(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="trend")
      tp3(g3, lk, gam(first1,0), iota(first1,0));
    else if(dynamics=="ricker")
      tp4(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="gompertz")
      tp5(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
  } else if(go_dims == "rowvec") {
    for(int t=0; t<(T-1); t++) {
      if(ytna(first1,t)==1) { // FIXME: this is not generic!
	      continue;
      }
      if(dynamics=="constant" || dynamics=="notrend") {
	      tp1(g3_t.slice(t), nrI, nrI1, N, I, I1, Ib, Ip, 
            gam(first1,t), om(first1,t));
      }
      else if(dynamics=="autoreg") {
	      tp2(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
      }
      else if(dynamics=="trend")
	      tp3(g3_t.slice(t), lk, gam(first1,t), iota(first1,t));
      else if(dynamics=="ricker")
	      tp4(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
      else if(dynamics=="gompertz")
	      tp5(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
    }
  }

  // loop over sites
  for(int i=0; i<M; i++) {
    ll_i=0.0;
    first_i = first[i]-1; // remember 0=1st location in C++
    last_i = last[i]-1;
    
    //Calculate g_star
    g_star.ones();
    if(last_i > first_i) {
      // loop over time periods in reverse order, up to second occasion
      for(int t=last_i; t>first_i; t--) {	      
        //Skip site i time t if all NA
        if(ytna(i,t)==1) {
	        continue;
	      }

        //Detection
        vec ysub = y.subcube(span(i), span(t), span());
        vec cp = distprob(keyfun, sig(i,t), scale, survey, 
                          db, w, a.row(i)) % u.col(i);
        double ycp = sum(ysub % log(cp));
        vec lkmyt_sub = lfac_kmyt.subcube(span(i), span(t), span());
             

        //vec p1 = y.subcube(span(i), span(t), span());
        //vec cp = distprob(keyfun, sig(i,t), scale, survey, 
                          //db, w, a.row(i));
        //vec p2 = cp % u.col(i);
        //vec p3 = lgy1(span(i), span(t), span());
        //vec p5 = lfac_kmyt.subcube(span(i), span(t), span());

        //double p2m3 = sum(p1 % log(p2));
        
        if(keyfun == "uniform"){
          g1_t.zeros();
          unsigned kint = sum(ysub);
          g1_t(kint) = exp(lfac_k(kint) - lkmyt_sub(kint) + ycp);

        } else {

          vec fin_sub = fin.subcube(span(i), span(t), span());
          vec kmyt_sub = kmyt.subcube(span(i), span(t), span());
          double cpJ = 1 - sum(cp);

          g1_t = lfac_k - lkmyt_sub + ycp + log(cpJ) * kmyt_sub;
          g1_t = exp(g1_t);
          g1_t = g1_t % fin_sub; //set to 0 when K impossible
        }

	      g1_t_star = g1_t % g_star;

        // computes transition probs for g3
        if(go_dims == "matrix") {
	        g3.zeros();
	        if(dynamics=="constant" || dynamics=="notrend") {
	          tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(i,t-1), om(i,t-1));
	        }
	        else if(dynamics=="autoreg") {
	          tp2(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	        }
	        else if(dynamics=="trend")
	          tp3(g3, lk, gam(i,t-1), iota(i,t-1));
	        else if(dynamics=="ricker")
	          tp4(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	        else if(dynamics=="gompertz")
	          tp5(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
        } 
        else if(go_dims == "rowvec") {
	        g3 = g3_t.slice(t-1);
        }
        int delta_it = delta(i,t);
        // matrix multiply transition probs over time gaps
        if(delta_it>1) {
	        g3_d = g3;
	      for(int d=1; d<delta_it; d++) {
	        g3_d = g3_d * g3;
	      }
	      g_star = g3_d * g1_t_star;
        } 
        else {
	        g_star = g3 * g1_t_star;
        }
      }
    }
    
    //Calculate g1
    int delta_i0 = delta(i,0);
    g1.zeros();


    vec ysub = y.subcube(span(i), span(first_i), span());
    vec cp = distprob(keyfun, sig(i,first_i), scale, survey, 
                          db, w, a.row(i)) % u.col(i);
    double ycp = sum(ysub % log(cp));
    vec lkmyt_sub = lfac_kmyt.subcube(span(i), span(first_i), span());

    //vec p1 = y.subcube(span(i), span(first_i), span());
    //vec cp = distprob(keyfun, sig(i,first_i), exp(scale), survey, 
                          //db, w, a.row(i));
    //vec p2 = cp % u.col(i);
    //vec p3 = lgy1(span(i), span(first_i), span());
    //double p2m3 = sum(p1 % log(p2) - p3);

    //double p2m3 = sum(p1 % log(p2));

    //vec p5 = lfac_kmyt.subcube(span(i), span(first_i), span());
    
    if(keyfun == "uniform"){
      g1.zeros();
      unsigned kint = sum(ysub);
      g1(kint) = exp(lfac_k(kint) - lkmyt_sub(kint) + ycp);

    } else {
      vec fin_sub = fin.subcube(span(i), span(first_i), span());
      vec kmyt_sub = kmyt.subcube(span(i), span(first_i), span());
      double cpJ = 1 - sum(cp);

      g1 = lfac_k - lkmyt_sub + ycp + log(cpJ) * kmyt_sub;
      g1 = exp(g1);
      g1 = g1 % fin_sub; //set to 0 when K impossible
    }

    if(delta_i0>1){
	    g1_star = g1_t % g_star;
    }

    //Calculate g2
    //Nest k loops inside if statements to save time (?)
    colvec g2 = zeros(lk);
    if(mixture=="P"){
      for(int k=0; k<lk; k++){
        g2(k) = Rf_dpois(k, lam(i), false);
      }
    }
    else if(mixture=="NB") {
      for(int k=0; k<lk; k++){
        g2(k) = dnbinom_mu(k, alpha, lam(i), false);
      }
    }
    else if(mixture=="ZIP") {
      for(int k=0; k<lk; k++){
        g2(k) = dzip(k, lam(i), psi);
      }
    }
    
    //Combine g1, g2, gstar
    if(delta_i0==1) {
      vec g12star = g1 % g2 % g_star; 
      ll_i += log(sum(g12star));     
    }
    else if(delta_i0>1) {
      g3_d = g3;
      for(int d=0; d<delta_i0; d++) {
        g3_d = g3_d * g3;
      }
      g_star = g3_d * g1_star;
      vec g2star = g2 % g_star;
      ll_i += log(sum(g2star));     
    }

    ll += ll_i;

  }

  return wrap(-ll);
}
