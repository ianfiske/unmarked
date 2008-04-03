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

//#define DEBUG


//typedef Matrix<double,Col,Concrete> matrix;

using namespace std;
using namespace scythe;
using std::vector;

class HiddenMarkovModel {
public:
    double operator() (const Matrix<> &parms);
    Matrix<> Hdet;
    Matrix<> Hphi;
    unsigned int nDMP, nDMP_un, nPhiP_un, nDYP, nDP, nSP, nPhiP,
	nP, K, nDCP, M, J, nY;
    bool yearly_det;
    vector< vector< vector <unsigned int> > > y_itj;

// Design matrix for each site, time, observation. dim:
// K x ncol(XDet)
    vector< vector< vector < Matrix<> > > > XDet_itj;
};


// This is the transpose of the phi matrix in my notes.  ie., it is column 
// stochastic
inline Matrix<> phiMatrix(const Matrix<>& phiParms)
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

// detParms is vector containing p1:p3, beta21, beta32, beta31.
inline Matrix<> detMatrix(const Matrix<>& detParms)
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

double  HiddenMarkovModel::operator() (const Matrix<>& parms)
{

#ifdef DEBUG
    cout << "parms:\n" << parms << endl;
#endif

    const Matrix<> detParms = Hdet * parms(0, 0, (nDMP - 1), 0);

#ifdef DEBUG
    cout << "detParms:\n" << detParms << endl;
#endif

    const Matrix<> alpha = detParms(0, 0, K - 1, 0);

#ifdef DEBUG
    cout << "alpha1" << endl;
#endif

    const Matrix<> logitBeta = detParms(K, 0, nDMP_un - 1, 0);
    const Matrix<> beta = plogis(logitBeta, 0, 1);
    Matrix<> detParmVec = alpha;

#ifdef DEBUG
    cout << "alpha:\n" << alpha << endl;
#endif

    Matrix<> detMatObs(K + 1, K + 1, false, 0);
    Matrix<> p(K, 1, false, 0);
    Matrix<> detVecObs;
    Matrix<> detVec;
    double negLogLike = 0;
    Matrix<> psiSite;

#ifdef DEBUG
    cout << "beta:\n" << beta << endl;
#endif
    
    if (nDCP > 0) {
	const Matrix<> b = parms(nDMP, 0, nDMP + nDCP - 1, 0);
	detParmVec = rbind(detParmVec, b);
    }

    if (yearly_det) {
	const Matrix<> gma = parms(nDMP + nDCP, 0, nDP - 1, 0);
	detParmVec = rbind(gma, detParmVec);
    }
  
    Matrix<> logPsi = parms(nDP, 0, nDP + K - 1, 0);
    Matrix<> Zero(0);
    logPsi = rbind(Zero, logPsi);
    Matrix<> psi = scythe::exp(logPsi);
    psi = psi / sum(psi);

#ifdef DEBUG
    cout << "psi:\n" << psi;
#endif

    const Matrix<> logitPhiParms = parms(nDP + K, 0, nP - 1, 0);
    const Matrix<> phiParms =  Hphi * plogis(logitPhiParms, 0, 1);
    const Matrix<> phi = phiMatrix(phiParms);

#ifdef DEBUG
    cout << "phi:\n" << phi << endl;
#endif

    for (unsigned int i = 0; i < M; ++i) {
	
	// Initialize psi for site i.
	psiSite = psi;

	for (unsigned int t = 0; t < nY; ++t) {
	    // initialize detVec to ones.
	    // this will be the product the observation's detVec's.
	    detVec = ones(K + 1, 1);
	    
	    for (unsigned int j = 0; j < J; ++j) {
		// compute stasis probabilities: p = plogis();
		p = plogis(XDet_itj[i][t][j] * detParmVec, 0, 1);

#ifdef DEBUG
		cout << "XDet:\n" << XDet_itj[i][t][j] 
		     << "\ndetParmVec:\n" << detParmVec
		     << "\np:\n" << p << endl;
#endif

		// get detection matrix for obs j
		detMatObs = detMatrix(rbind(p, beta));

		// pull out the proper column, according to y
		detVecObs = detMatObs(_,y_itj[i][t][j]);
#ifdef DEBUG
		cout << "detMatObs:\n" << detMatObs
		     << "\ny:" << y_itj[i][t][j]
		     << "\ndetVecObs:\n" << detVecObs << endl;
#endif

		// update yearly detection vector;
		detVec = detVec % detVecObs;
	    }

	    // mult by appropriate detMat and phi 
	    psiSite = psiSite % detVec;

	    // if last time, then instead of multiplying by phi^t, just
	    // sum up the vector.  this is the likelihood for site i.
	    if (t < (nY - 1))
		psiSite = phi * psiSite;
	    else
		negLogLike += -log(sum(psiSite));

#ifdef DEBUG
	    cout << "detVec:\n" << detVec
		 << "psiSite:\n" << psiSite << endl;
#endif

	}
    }

#ifdef DEBUG
    cout << "negll: " << negLogLike << endl;
#endif

    return negLogLike;
}
    
extern "C" {
    
    void findMLE (int *y_itj,
		  double *XDet_itjk,
		  const int *ncXDet,
		  const int *nDMP, 
		  const int *nDCP, 
		  const int *nDP, 
		  const int *nDYP,
		  const int *nSP, 
		  const int *nPhiP, 
		  const int *nP,
		  const int *nDMP_un, 
		  const int *nPhiP_un,
		  const double *Hdet_in,
		  const double *Hphi_in, 
		  const int *K,
		  const int *yearly_det, 
		  const int *M,
		  const int *J,
		  const int *nY,
		  double *MLE)
    {
	#ifdef DEBUG
	cout << "entered findMLE" << endl;
	#endif

	// set random number generator
	mersenne myrng;
	
	// instantiate model
	HiddenMarkovModel hmm;
	
	#ifdef DEBUG
	cout << "model hmm had been instantiated" << endl;
	#endif

	for (int i = 0; i < *M; ++i) {
	    hmm.y_itj.push_back(vector< vector< unsigned int > > () );
	    hmm.XDet_itj.push_back(vector< vector< Matrix<> > > () );
	    for (int t = 0; t < *nY; ++t) {
		hmm.y_itj[i].push_back(vector< unsigned int > () );
		hmm.XDet_itj[i].push_back(vector< Matrix<> > () );
		for (int j = 0; j < *J; ++j) {
		    hmm.y_itj[i][t].push_back(*(y_itj++));
		    // note that this is the tranpose of XDet
		    hmm.XDet_itj[i][t].push_back(Matrix<> (*ncXDet, *K, XDet_itjk));
		    XDet_itjk += (*K) * (*ncXDet);  // this line will change XDet changes.
		    hmm.XDet_itj[i][t][j] = scythe::t(hmm.XDet_itj[i][t][j]);
		}
	    }
	}

	#ifdef DEBUG
	cout << "Model data has been input." << endl;
	#endif

	hmm.Hdet = Matrix<> (*nDMP_un, *nDMP, Hdet_in);
	hmm.Hphi = Matrix<> (*nPhiP_un, *nPhiP, Hphi_in);
	hmm.nDMP = *nDMP;
	hmm.nDCP = *nDCP;
	hmm.nDP = *nDP;
	hmm.nDYP = *nDYP;
	hmm.nSP = *nSP;
	hmm.nPhiP = *nPhiP;
	hmm.nP = *nP;
	hmm.nDMP_un = *nDMP_un;
	hmm.nPhiP_un = *nPhiP_un;
	hmm.K = *K;
	hmm.yearly_det = *yearly_det;
	hmm.M = *M;
	hmm.J = *J;
	hmm.nY = *nY;
	
	#ifdef DEBUG
	cout << "Hdet:\n" << hmm.Hdet << endl;
	#endif
	    

	
	Matrix<> theta(*nP, 1, true, 0);

	Matrix<> theta_MLE;
	try{
	    theta_MLE = scythe::BFGS(hmm, theta, myrng, 1e5, 1e-5, false);
	} catch (scythe_convergence_error err) {
	    cout << "BFGS did not converge." << endl;	    
	}

#ifdef DEBUG
	cout << "MLE:\n" << theta_MLE << endl;
#endif
	for (int i = 0; i < *nP; ++i)
	    MLE[i] = theta_MLE(i);

    }
}
