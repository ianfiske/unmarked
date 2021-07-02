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

//Integrate with trapezoidal rule
double trap_rule(IntFunc &f, double a, double b){
  
  int n = 100;
  double h = (b - a) / n;
  
  double int_sum = 0;
  for(int i=1; i<n; i++){
    int_sum += f(a + i*h);
  }

  return( h/2 * (f(a) + 2*int_sum + f(b)) );
}
