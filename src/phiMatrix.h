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


/* Create a class for phi matrices (functors). Then, during model instantiation,
   set the phi matrix to the appropriate functor for that model based on K.
   Maybe function pointers are better?
*/

// This is the transpose of the phi matrix in my notes.  ie., it is column 
// stochastic
using namespace std;
using namespace scythe;
using std::vector;

class phiMatrix {
public:
    virtual Matrix<> operator() (const Matrix<>& phiParms) =0;
};

class phiMatrix4 : public phiMatrix {
public:
    Matrix<> operator() (const Matrix<>& phiParms)
    {
	double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParms(2), s3 = phiParms(3),
	    g01 = phiParms(4), g02 = phiParms(5), g12 = phiParms(6), 
	    r10 = phiParms(7), r21 = phiParms(8), r20 = phiParms(9),
	    r32 = phiParms(10), r31 = phiParms(11);
	double vals[16] = {s0, (1-s0)*g01, (1-s0)*(1-g01)*g02, (1-s0)*(1-g01)*(1-g02),
			   (1-s1)*(1-g12)*r10, s1, (1-s1)*g12, (1-s1)*(1-g12)*(1-r21),
			   (1-s2)*(1-r21)*r20, (1-s2)*r21, s2, (1-s2)*(1-r21)*(1-r20),
			   (1-s3)*(1-r32)*(1-r31), (1-s3)*(1-r32)*r31, (1-s3)*r32, s3};
    
	Matrix<> P(4, 4, vals);
	return P;
    }
};

class phiMatrix3 : public phiMatrix {
public:
    Matrix<> operator() (const Matrix<>& phiParms)
    {
	double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParms(2),
	    g0 = phiParms(3), g1 = phiParms(4), r2 = phiParms(5);

	double vals[9] = {s0,        (1-s0)*g0,     (1-s0)*(1-g0),
			  (1-s1)*(1-g1), s1,        (1-s1)*g1,
			  (1-s2)*(1-r2), (1-s2)*r2,  s2};

	Matrix<> P(3, 3, vals);

	return P;
    }
};

