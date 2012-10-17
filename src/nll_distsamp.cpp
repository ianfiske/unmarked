#include "nll_distsamp.h"
#include "detfuns.h"


SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP scale_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP survey_, SEXP reltol_ ) {

  Rcpp::IntegerMatrix y(y_);
  Rcpp::NumericVector lam(lam_);
  Rcpp::NumericVector sig(sig_);
  double scale = Rcpp::as<double>(scale_);
  Rcpp::NumericMatrix a(a_);
  Rcpp::NumericMatrix u(u_);
  Rcpp::NumericVector w(w_);
  Rcpp::NumericVector db(db_);
  std::string keyfun = Rcpp::as<std::string>(keyfun_);
  std::string survey = Rcpp::as<std::string>(survey_);

  bool verbose=false;

  int R = y.nrow();   //y.n_rows;
  int J = y.ncol();   // y.n_cols;

  //  double mu = 0.0;
  double ll = 0.0;
  double lnmin = log(DOUBLE_XMIN);

  //  void *ex;
  double *ex;
  ex = (double *) R_alloc(2, sizeof(double));
  // Integration settings given to Rdqags
  double lower = 0.0, upper = 0.0, epsrel = 0.0, epsabs = 0.0;
  int limit = 100, lenw = 400, last=0, iwork = 100;
  double work = 0.0, result = 0.0, abserr = 0.0;
  int neval = 0, ier=0;
  double cp = 0.0;

  double f0 = 0.0;

  for(int i=0; i<R; i++) {
    if((survey=="line") & (keyfun=="halfnorm"))
      f0 = Rf_dnorm4(0.0, 0.0, sig[i], false);
    if((survey=="line") & (keyfun=="exp"))
      f0 = Rf_dexp(0.0, 1/sig[i], false);
    for(int j=0; j<J; j++) {
      cp = 0.0;
      ex[0] = sig[i];
      ex[1] = scale;
      // Integration settings given to Rdqags
      lower = db[j];
      upper = db[j+1];
      epsrel = Rcpp::as<double>(reltol_);
      epsabs = epsrel;
      limit = 100;
      lenw = 400;
      last = 0;
      iwork = 100;
      work = 400.0;
      result = DOUBLE_XMIN;
      abserr = 0.0;
      neval = 0;
      ier = 0;
      if(keyfun=="uniform") {
	cp = u(i,j);
      } else {
	if(survey=="point") {
	  if(keyfun=="halfnorm") {
	    result = sig[i]*sig[i]*(1-exp(-upper*upper/(2*sig[i]*sig[i])))-
	      sig[i]*sig[i]*(1-exp(-lower*lower/(2*sig[i]*sig[i])));
	    // Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    // 	   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    // 	   &work);
	  } else if(keyfun=="exp") {
	    Rdqags(grexp, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  } else if(keyfun=="hazard") {
	    Rdqags(grhaz, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  }
	  if(ier > 0 && verbose) {
	    Rf_warning("The integration was not successful.");
	  }
	  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	} else if(survey=="line") {
	  if(keyfun=="halfnorm") {
	    result = (Rf_pnorm5(upper, 0.0, sig[i], true, false) -
		      Rf_pnorm5(lower, 0.0, sig[i], true, false)) / f0;
	    // Rdqags(gxhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    //     &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    //	   &work);
	  } else if(keyfun=="exp") {
	    result = (Rf_pexp(upper, 1/sig[i], true, false) -
		      Rf_pexp(lower, 1/sig[i], true, false)) / f0;
	    // Rdqags(gxexp, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    // 	   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    // 	   &work);
	  } else if(keyfun=="hazard") {
	    Rdqags(gxhaz, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  }
	  if(ier > 0 && verbose) {
	    Rf_warning("Warning: the integration was not successful.");
	  }
	  cp = result / w[j] * u(i,j);
	}
      }
      ll += std::max(Rf_dpois(y(i,j), lam[i]*cp, true), lnmin);
    }
  }
  return Rcpp::wrap(-ll);
}
