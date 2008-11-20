#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "lapack.h"
#include "ide.h"
#include "smath.h"
#include "mersenne.h"
#include "optimize.h"
#include <iostream>
#include <R.h>
#include <vector>

using namespace std;
using namespace scythe;
using std::vector;

class detMatrix {
public:
  //    virtual Matrix<> operator() (const Matrix<>& detParms) =0;
    virtual Matrix<> operator() (const Matrix<>& detParms, const int& y) =0;
};
/*
// detParms is vector containing p1:p3, beta21, beta32, beta31.
class detMatrix4: public detMatrix {
public:
    Matrix<> operator() (const Matrix<>& detParms)
    {
	double p1 = detParms(0), p2 = detParms(1), p3 = detParms(2),
	    beta21 = detParms(3), beta32 = detParms(4), beta31 = detParms(5);

	double vals[16] = {1, 1-p1, (1-beta21)*(1-p2),  (1-p3)*(1-beta32)*(1-beta31),
			   0,  p1,    beta21*(1-p2),      beta31*(1-beta32)*(1-p3),
			   0,  0,      p2,                 beta32*(1-p3),
			   0,  0,      0,                      p3};

	Matrix<> D(4, 4, vals);
	return D;
    }
};


class detMatrix3: public detMatrix {
public:
    Matrix<> operator() (const Matrix<>& detParms)
    {
	double p1 = detParms(0), p2 = detParms(1), beta2 = detParms(2);
	double vals[9] = {1, 1-p1, (1-p2)*(1-beta2),
			  0, p1,     (1-p2)*beta2,
			  0, 0,        p2};
	Matrix<> D(3, 3, vals);
	return D;
    }
};

class detMatrix2: public detMatrix {
public:
    Matrix<> operator() (const Matrix<>& detParms)
    {
	double p = detParms(0);

	double vals[4] = {1,   1-p,
			  0,    p};

	Matrix<double,Col,View> D(2, 2, vals);


	return D;
    }
};
*/

class detLogit2: public detMatrix {

 public:
  Matrix<> operator() (const Matrix<>& detParms, const int& y) {

    double p = detParms(0);
    Matrix<> dv(2, 1, false);
    if(y == 0) {
      double vals[2] = {1, 1/(1 + exp(p))};
      dv = Matrix<>(2,1,vals);
    } else {
      double vals[2] = {0, exp(p)/(1 + exp(p))};
      dv = Matrix<>(2,1,vals);
    }
    return dv;
  }

};

class detLogit3: public detMatrix {
public:
  Matrix<> operator() (const Matrix<>& detParms, const int& y)
    {

      double p0 = detParms(0), p1 = detParms(1), p2 = detParms(2);

      Matrix<> dv(3,1,false);
      switch (y) {
      case 0:
	{
	  double vals[3] = {1, 1/(1 + exp(p0)),
			    1/(1 + exp(p1) + exp(p2))};
	  dv = Matrix<>(3, 1, vals);
	}
	break;
      case 1:
	{
	  double vals[3] = {0, exp(p0)/(1 + exp(p0)),
			    exp(p1)/(1 + exp(p1) + exp(p2))};
	  dv = Matrix<>(3, 1, vals);
	}
	break;
      case 2:
	{
	  double vals[3] = {0, 0, exp(p2)/(1 + exp(p1) + exp(p2))};
	  dv = Matrix<>(3, 1, vals);
	}
	break;
      }
			      
      return dv;
    }
};


// detParms is vector containing p1:p3, beta21, beta32, beta31.
class detLogit4: public detMatrix {
 public:
  Matrix<> operator() (const Matrix<>& detParms, const int& y)
    {

      double p1 = detParms(0), p2 = detParms(1), p3 = detParms(2),
	p4 = detParms(3), p5 = detParms(4), p6 = detParms(5);

      Matrix<> dv(4,1,false);
      switch (y) {
      case 0:
	{
	  double vals[4] = {1, exp(p1)/(1 + exp(p1)),
			    exp(p2)/(1 + exp(p2) + exp(p3)),
			    exp(p4)/(1 + exp(p4) + exp(p5) + exp(p6))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 1:
	{
	  double vals[4] = {0, 1/(1 + exp(p1)),
			    exp(p3)/(1 + exp(p2) + exp(p3)),
			    exp(p5)/(1 + exp(p4) + exp(p5) + exp(p6))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 2:
	{
	  double vals[4] = {0, 0, 1/(1 + exp(p2) + exp(p3)),
			    exp(p6)/(1 + exp(p4) + exp(p5) + exp(p6))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 3:
	{
	  double vals[4] = {0, 0, 0, 
			    1/(1 + exp(p4) + exp(p5) + exp(p6))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      }
			      

      //	Matrix<> D(4, 4, vals);
      //	return scythe::t(D);
      return dv;
    }
};

class detLogit4ar: public detMatrix {
 public:
  Matrix<> operator() (const Matrix<>& detParms, const int& y)
    {

      double r1 = detParms(0), r2 = detParms(1), r3 = detParms(2);

      Matrix<> dv(4,1,false);
      switch (y) {
      case 0:
	{
	  double vals[4] = {1, exp(r1)/(1 + exp(r1)),
			    exp(2*r2)/(1 + exp(2*r2) + exp(r2)),
			    exp(3*r3)/(1 + exp(3*r3) + exp(2*r3) + exp(r3))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 1:
	{
	  double vals[4] = {0, 1/(1 + exp(r1)),
			    exp(r2)/(1 + exp(2*r2) + exp(r2)),
			    exp(2*r3)/(1 + exp(3*r3) + exp(2*r3) + exp(r3))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 2:
	{
	  double vals[4] = {0, 0, 1/(1 + exp(2*r2) + exp(r2)),
			    exp(r3)/(1 + exp(3*r3) + exp(2*r3) + exp(r3))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      case 3:
	{
	  double vals[4] = {0, 0, 0, 
			    1/(1 + exp(3*r3) + exp(2*r3) + exp(r3))};
	  dv = Matrix<>(4, 1, vals);
	}
	break;
      }
			      

      //	Matrix<> D(4, 4, vals);
      //	return scythe::t(D);
      return dv;
    }
};
