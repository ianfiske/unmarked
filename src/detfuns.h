#ifndef _UNMARKED_DETFUNS_H
#define _UNMARKED_DETFUNS_H

#include <Rcpp.h>

using namespace Rcpp;

//Base class for functions to integrate
class IntFunc {
  public:
    IntFunc() {}
    
    virtual double operator()(const double& x) const {
      return(x);
    }
};

//Negative exponential function
class DetExp: public IntFunc {
  private:
    double rate;
    int point;
  public:
    DetExp(double rate_, int point_) :
      rate(rate_), point(point_) {}

    double operator()(const double& x) const {
      double pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( std::exp( -x/rate) * pd_adjust);
    }
};

//Hazard function
class DetHaz: public IntFunc {
  private:
    double shape;
    double scale;
    int point;
  public:
    DetHaz(double shape_, double scale_, int point_) : 
      shape(shape_), scale(scale_), point(point_) {}

    double operator()(const double& x) const {
      double pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( (1-std::exp(-1 * pow(x/shape, -scale))) * pd_adjust);
    }
};

double trap_rule(IntFunc &f, double a, double b) ;

#endif
