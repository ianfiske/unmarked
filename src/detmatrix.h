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
    virtual Matrix<> operator() (const Matrix<>& detParms) =0;
};

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
