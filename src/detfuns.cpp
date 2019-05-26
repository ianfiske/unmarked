#include "detfuns.h"

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
