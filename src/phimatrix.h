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
//#include <R.h>
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

// Enumeration value defitions for phi matrix string choices
enum phiString { evNotDefined,
		 evMacKenzie,
		 evThreeState,
		 evFourState,
		 evFourState4,
		 evFourStateAR,
		 evFourStateRandom,
		 evCumLogit,
		 evLogit4,
		 evLogit4ar,
		 evLogit3,
		 evLogit2,
		 evEnd };

// Map to associate the strings with the enum values
static std::map<std::string, phiString> phi_mapStringValues;


class phiMatrix {
public:
    virtual Matrix<> operator() (const Matrix<>& phiParms) =0;
};

/* class phiMatrix4 : public phiMatrix { */
/* public: */
/*     Matrix<> operator() (const Matrix<>& phiParms) */
/*     { */
/* 	const double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParms(2), s3 = phiParms(3), */
/* 	    g01 = phiParms(4), g02 = phiParms(5), g12 = phiParms(6), */
/* 	    r10 = phiParms(7), r21 = phiParms(8), r20 = phiParms(9), */
/* 	    r32 = phiParms(10), r31 = phiParms(11); */
/* 	const double vals[16] = {s0, (1-s0)*g01, (1-s0)*(1-g01)*g02, (1-s0)*(1-g01)*(1-g02), */
/* 			   (1-s1)*(1-g12)*r10, s1, (1-s1)*g12, (1-s1)*(1-g12)*(1-r10), */
/* 			   (1-s2)*(1-r21)*r20, (1-s2)*r21, s2, (1-s2)*(1-r21)*(1-r20), */
/* 			   (1-s3)*(1-r32)*(1-r31), (1-s3)*(1-r32)*r31, (1-s3)*r32, s3}; */

/* 	Matrix<> P(4, 4, vals); */
/* 	return P; */
/*     } */
/* }; */

/* class phiMatrix3 : public phiMatrix { */
/* public: */
/*     Matrix<> operator() (const Matrix<>& phiParms) */
/*     { */
/* 	double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParms(2), */
/* 	    g0 = phiParms(3), g1 = phiParms(4), r2 = phiParms(5); */

/* 	double vals[9] = {s0,        (1-s0)*g0,     (1-s0)*(1-g0), */
/* 			  (1-s1)*(1-g1), s1,        (1-s1)*g1, */
/* 			  (1-s2)*(1-r2), (1-s2)*r2,  s2}; */

/* 	Matrix<> P(3, 3, vals); */

/* 	return P; */
/*     } */
/* }; */

/* class phiMatrix2 : public phiMatrix { */
/* public: */
/*     Matrix<> operator() (const Matrix<>& phiParms) */
/*     { */
/* 	double s0 = phiParms(0), s1 = phiParms(1); */

/* 	double vals[4] = {s0,    1-s0, */
/* 			  1-s1,   s1}; */

/* 	Matrix<> P(2, 2, vals); */

/* 	return P; */
/*     } */
/* }; */


/* class phiMatrix4p4 : public phiMatrix { */
/*  public: */
/*   Matrix<> operator() (const Matrix<>& phiParmsIn) */
/*     { */
/*       Matrix<> phiParms = plogis(phiParmsIn, 0, 1); */

/*       double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParm, */
/* 	g = phiParms(2), r = phiParms(3); */

/*       double vals[16] = {s0,    g,  g*g,   g*g*g, */
/* 			 r,     s1,   g,     g*g, */
/* 			 r*r,    r,  s2,      g, */
/* 			 r*r*r, r*r, r,      s3}; */

/*       const Matrix<> P(4, 4, vals); */

/*       const Matrix<> Pcolsums = ones(1,4) * P; */

/*       for(int i = 0; i < 4; ++i) */
/* 	P(_,i) = P(_,i)/Pcolsums(i); */

/*       return(P); */
/*     } */
/* }; */

class phiMatrix4ar : public phiMatrix {
 public:
  Matrix<> operator() (const Matrix<>& phiParmsIn)
    {

      Matrix<> phiParms = plogis(phiParmsIn, 0, 1);

      double s0 = phiParms(0), s1 = phiParms(1), s2 = phiParms(2),
	s3 = phiParms(3), g0 = phiParms(4), g1 = phiParms(5),
	g2 = phiParms(6), r1 = phiParms(7), r2 = phiParms(8),
	r3 = phiParms(9);

      double vals[16] = {s0,    g0,    g0*g0,   g0*g0*g0,
			 r1,     s1,     g1,     g1*g1,
			 r2*r2,    r2,    s2,      g2,
			 r3*r3*r3, r3*r3, r3,      s3};

      const Matrix<> P(4, 4, vals);

      const Matrix<> Pcolsums = ones(1,4) * P;

      for(int i = 0; i < 4; ++i)
	P(_,i) = P(_,i)/Pcolsums(i);

      return(P);
    }
};

/*class cumLogitPhiMatrix : public phiMatrix {
public:
    Matrix<> operator() (const Matrix<>& phiParms)
    {
	// There are K alphas (none for last state)
	// phiParm(0) is baseline alpha
	// phiParm(j) is sqrt(alpha(j) - alpha(j-1)) so that
	// alphas are increasing in j
//	Matrix<> alpha = phiParms(0, 0, K - 1, 0);
	double alpha[K];
	double beta = phiParms(K,0);

	alpha[0] = phiParms(0);
	Matrix<> alphaDiffs = phiParms(1, 0, K - 1, 0);
	alphaDiffs = alphaDiffs^2;

	for(int j = 1; j <= K - 1; ++j) {
	    alpha[j] = alpha[j - 1] + alphaDiffs(j - 1);
	}

#ifdef DEBUG
	cout << "phiParms:\n" << phiParms
	     << "alphaDiffs:\n" << alphaDiffs
	     << "\nalpha: \n" << Matrix<> (K, 1, alpha)
	     << "\nbeta: " << beta << endl;
#endif
	Matrix<> P(K + 1, K + 1, true, 0);

	// assign values columnwise
	// element (j,i) is Prob(Xt = j| Xt-1 = i)
	for(int i = 0; i <= K; i++) {
	    P(0,i) = plogis(alpha[0] + beta*i, 0, 1);
	    for(int j = 1; j < K; j++) {
		P(j,i) = plogis(alpha[j] + beta*i, 0, 1) -
		    plogis(alpha[j-1] + beta*i, 0, 1);
	    }
	    P(K,i) = 1 - sum(P(_,i));
	}
#ifdef DEBUG
	cout << "P:\n" << P << endl;
#endif
	return(P);
    }

    int K;
};
*/

class phiLogit2 : public phiMatrix {
 public:
  Matrix<> operator() (const Matrix<>& phiParms) {

    const double p0 = phiParms(0), p1 = phiParms(1);

    const double vals[4] = {1/(1 + exp(p0)), exp(p0)/(1 + exp(p0)),
			    1/(1 + exp(p1)), exp(p1)/(1 + exp(p1))};

    Matrix<> P(2, 2, vals);

    return(P);

  }

};

class phiLogit3 : public phiMatrix {
 public:
  Matrix<> operator() (const Matrix<>& phiParms) {

    const double p0 = phiParms(0), p1 = phiParms(1), p2 = phiParms(2),
      p3 = phiParms(3), p4 = phiParms(4), p5 = phiParms(5);


    const double vals[9] = {1/(1 + exp(p0) + exp(p1)),
			     exp(p0)/(1 + exp(p0) + exp(p1)),
			     exp(p1)/(1 + exp(p0) + exp(p1)),
			     1/(1 + exp(p2) + exp(p3)),
			     exp(p2)/(1 + exp(p2) + exp(p3)),
			     exp(p3)/(1 + exp(p2) + exp(p3)),
			     1/(1 + exp(p4) + exp(p5)),
			     exp(p4)/(1 + exp(p4) + exp(p5)),
			     exp(p5)/(1 + exp(p4) + exp(p5))};


    const Matrix<> P(3, 3, vals);

    return(P);

  }
};

class phiLogit4 : public phiMatrix {
 public:
  Matrix<> operator() (const Matrix<>& phiParms) {

    const double p0 = phiParms(0), p1 = phiParms(1), p2 = phiParms(2),
      p3 = phiParms(3), p4 = phiParms(4), p5 = phiParms(5),
      p6 = phiParms(6), p7 = phiParms(7), p8 = phiParms(8),
      p9 = phiParms(9), p10 = phiParms(10), p11 = phiParms(11);

    const double vals[16] = {1/(1 + exp(p0) + exp(p1) + exp(p2)),
			     exp(p0)/(1 + exp(p0) + exp(p1) + exp(p2)),
			     exp(p1)/(1 + exp(p0) + exp(p1) + exp(p2)),
			     exp(p2)/(1 + exp(p0) + exp(p1) + exp(p2)),
			     1/(1 + exp(p3) + exp(p4) + exp(p5)),
			     exp(p3)/(1 + exp(p3) + exp(p4) + exp(p5)),
			     exp(p4)/(1 + exp(p3) + exp(p4) + exp(p5)),
			     exp(p5)/(1 + exp(p3) + exp(p4) + exp(p5)),
			     1/(1 + exp(p6) + exp(p7) + exp(p8)),
			     exp(p6)/(1 + exp(p6) + exp(p7) + exp(p8)),
			     exp(p7)/(1 + exp(p6) + exp(p7) + exp(p8)),
			     exp(p8)/(1 + exp(p6) + exp(p7) + exp(p8)),
			     1/(1 + exp(p9) + exp(p10) + exp(p11)),
			     exp(p9)/(1 + exp(p9) + exp(p10) + exp(p11)),
			     exp(p10)/(1 + exp(p9) + exp(p10) + exp(p11)),
			     exp(p11)/(1 + exp(p9) + exp(p10) + exp(p11))};


    const Matrix<> P(4, 4, vals);

    return(P);

  }
};


class phiLogit4ar : public phiMatrix {
 public:
  Matrix<> operator()(const Matrix<>& phiParms) {
    const Matrix<> expPhiParms = exp(phiParms);
    const double eb0 = expPhiParms(0), eb1 = expPhiParms(1), eb2 = expPhiParms(2),
      ea1 = expPhiParms(3), ea2 = expPhiParms(4), ea3 = expPhiParms(5);

    const double d0 = 1 + eb0 + pow(eb0,2) + pow(eb0,3);
    const double d1 = ea1 + 1 + eb1 + pow(eb1,2);
    const double d2 = pow(ea2,2) + ea2 + 1 + eb2;
    const double d3 = pow(ea3,3) + pow(ea3,2) + ea3 + 1;

    const double vals[16] = {1/d0, eb0/d0, pow(eb0,2)/d0, pow(eb0,3)/d0,
			     ea1/d1, 1/d1, eb1/d1, pow(eb1,2)/d1,
			     pow(ea2,2)/d2, ea2/d2, 1/d2, eb2/d2,
			     pow(ea3,3)/d3, pow(ea3,2)/d3, ea3/d3, 1/d3};

    const Matrix<> P(4, 4, vals);

    return(P);

  }
};

class phi3Absorb : public phiMatrix {
 public:
  Matrix<> operator()(const Matrix<>& phiParms) {
    const double gamma = plogis(phiParms(0), 0, 1),
      phi = plogis(phiParms(1), 0, 1);

    const double vals[9] = { 1 - gamma, gamma, 0,
			     0, phi, 1 - phi,
			     0, 0, 1};

    const Matrix<> P(3, 3, vals);
    return(P);
  }
};

class phi4random : public phiMatrix {
 public:
  Matrix<> operator() (const Matrix<>& phiParms) {
    const double p1 = exp(phiParms(0)),
      p2 = exp(phiParms(1)), p3 = exp(phiParms(2));
    const double d = 1 + p1 + p2 + p3;

    const double vals[16] = {1/d, p1/d, p2/d, p3/d,
			     1/d, p1/d, p2/d, p3/d,
			     1/d, p1/d, p2/d, p3/d,
			     1/d, p1/d, p2/d, p3/d};

    const Matrix<> P(4, 4, vals);
    return(P);
  }
};
