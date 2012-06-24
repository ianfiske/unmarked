#include "nll_distsamp.h"
#include "detfuns.h"


SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP survey_, SEXP reltol_ ) {

  Rcpp::IntegerMatrix y(y_);
  Rcpp::NumericVector lam(lam_);
  Rcpp::NumericVector sig(sig_);
  Rcpp::NumericMatrix a(a_);
  Rcpp::NumericMatrix u(u_);
  Rcpp::NumericVector w(w_);
  Rcpp::NumericVector db(db_);
  std::string keyfun = Rcpp::as<std::string>(keyfun_);
  std::string survey = Rcpp::as<std::string>(survey_);

  int R = y.nrow();   //y.n_rows;
  int J = y.ncol();   // y.n_cols;

  // Integration settings given to Rdqags
  // TODO: Allow for user-defined settings
  double epsabs = 0.01; // should be specific to measurement units
  double epsrel = Rcpp::as<double>(reltol_);
  int limit = 100;
  int lenw = 400;
  int last = 0;
  int iwork = 100;
  double work = 400.0;

  double mu = 0.0;
  double ll = 0.0;

  for(int i=0; i<R; i++) {
    for(int j=0; j<J; j++) {
      double cp = 0.0;
      if(keyfun=="uniform") {
	cp = u(i,j);
      } else {
	void *ex;
	ex = &sig(i);
	double lower = db[j];
	double upper = db[j+1];
	double result = 0.0;
	double abserr = 0.0;
	int neval = 0;
	int ier = 0;
	if(survey=="point") {
	  if(keyfun=="halfnorm") {
	    Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  } else if(keyfun=="exp") {
	    Rdqags(grexp, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  }
	  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	} else if(survey=="line") {
	  if(keyfun=="halfnorm") {
	    Rdqags(gxhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  } else if(keyfun=="exp") {
	    Rdqags(gxexp, ex, &lower, &upper, &epsabs, &epsrel, &result,
		   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
		   &work);
	  }
	  cp = result / w(j) * u(i,j);
	}
      }
      /* add error checking/handling here */
      //      mu = lam(i)*cp + DOUBLE_XMIN;
      ll += Rf_dpois(y(i,j), lam(i)*cp, true);
    }
  }
  return Rcpp::wrap(-ll);
}
