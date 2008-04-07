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
#include "phimatrix.h"
#include "detmatrix.h"

//#define DEBUG


//typedef Matrix<double,Col,Concrete> matrix;

using namespace std;
using namespace scythe;
using std::vector;

// base class for the library
class HiddenMarkovModel {
public:
    double operator() (const Matrix<> &parms);
    void fitmodel (const Matrix<> &inits);
    void loadData (int *y_itj_in, double *XDet_itjk, const int *ncXDet);
    void setMatrices (phiMatrix *phiMat, detMatrix *detMat);

    // constructor
    HiddenMarkovModel( const int *nDMP, 
		       const int *nDCP, 
		       const int *nDP, 
		       const int *nDYP,
		       const int *nSP, 
		       const int *nPhiP, 
		       const int *nP,
		       const int *nDMP_un, 
		       const int *nPhiP_un,
		       const int *K,
		       const int *M,
		       const int *J,
		       const int *nY,
		       const Matrix<> Hdet,
		       const Matrix<> Hphi,
		       const int *yearly_det):
	Hdet(Hdet), Hphi(Hphi),
	nDMP(*nDMP), nDCP(*nDCP), nDP(*nDP), nDYP(*nDYP), nSP(*nSP),
	nPhiP(*nPhiP), nP(*nP), nDMP_un(*nDMP_un), nPhiP_un(*nPhiP_un), 
	K(*K), M(*M), J(*J), nY(*nY),
	yearly_det(*yearly_det) { }

    Matrix<> Hdet;
    Matrix<> Hphi;
    unsigned int nDMP, nDCP, nDP, nDYP, nSP, nPhiP, nP, nDMP_un, nPhiP_un,
	K, M, J, nY;
    bool yearly_det;
    vector< vector< vector <unsigned int> > > y_itj;

// Design matrix for each site, time, observation. dim:
// K x ncol(XDet)
    vector< vector< vector < Matrix<> > > > XDet_itj;

    Matrix<> mle; // will store the MLE after fitting
    Matrix<> hessian;
    Matrix<> phi;

    phiMatrix *phiMatPtr;  // phi matrix functor pointer
    detMatrix *detMatPtr;  // det matrix functor pointer
};

void HiddenMarkovModel::loadData(int *y_itj_in,
				 double *XDet_itjk,
				 const int *ncXDet)
{
    for (unsigned int i = 0; i < M; ++i) {
	y_itj.push_back(vector< vector< unsigned int > > () );
	XDet_itj.push_back(vector< vector< Matrix<> > > () );
	for (unsigned int t = 0; t < nY; ++t) {
	    y_itj[i].push_back(vector< unsigned int > () );
	    XDet_itj[i].push_back(vector< Matrix<> > () );
	    for (unsigned int j = 0; j < J; ++j) {
		y_itj[i][t].push_back(*(y_itj_in++));
		    
		// note that this is the tranpose of XDet
		XDet_itj[i][t].push_back(Matrix<> (*ncXDet, K, XDet_itjk));
		    
		XDet_itjk += (K) * (*ncXDet);  // this line will change XDet changes.
		XDet_itj[i][t][j] = scythe::t(XDet_itj[i][t][j]);
	    }
	}
    }
}

void HiddenMarkovModel::setMatrices (phiMatrix *phiMatPtr, 
				     detMatrix *detMatPtr)
{
    this->phiMatPtr = phiMatPtr;
    this->detMatPtr = detMatPtr;
}

double HiddenMarkovModel::operator() (const Matrix<>& parms)
{

#ifdef DEBUG
    cout << parms << endl;
#endif

    const Matrix<> detParms = Hdet * parms(0, 0, (nDMP - 1), 0);
    const Matrix<> alpha = detParms(0, 0, K - 1, 0);
#ifdef DEBUG
    cout << "detParms:\n" << detParms << endl;
#endif
    const Matrix<> logitBeta = detParms(K, 0, nDMP_un - 1, 0);
    const Matrix<> beta = plogis(logitBeta, 0, 1);
    Matrix<> detParmVec = alpha;

    Matrix<> detMatObs(K + 1, K + 1, false, 0);
    Matrix<> p(K, 1, false, 0);
    Matrix<> detVecObs;
    Matrix<> detVec;
    double negLogLike = 0;
    Matrix<> psiSite;
#ifdef DEBUG
    cout << "alpha:\n" << alpha << endl;
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
    cout << "psi\n" << psi << endl;
#endif

    const Matrix<> logitPhiParms = parms(nDP + K, 0, nP - 1, 0);
    const Matrix<> phiParms =  Hphi * plogis(logitPhiParms, 0, 1);
    const Matrix<> phi = (*phiMatPtr)(phiParms);

#ifdef DEBUG
    cout << "phiParms\n" << phiParms
	 <<"\nphi\n" << phi << endl;
#endif

    for (vector<vector<vector< unsigned int > > >::size_type i = 0; i < M; ++i) {
	
	// Initialize psi for site i.
	psiSite = psi;

	for (unsigned int t = 0; t < nY; ++t) {
	    // initialize detVec to ones.
	    // this will be the product the observation's detVec's.
	    detVec = ones(K + 1, 1);
	    
	    for (unsigned int j = 0; j < J; ++j) {
		// 99 = missing value
		if (y_itj[i][t][j] != 99) {
		    // compute stasis probabilities: p = plogis();
		    p = plogis(XDet_itj[i][t][j] * detParmVec, 0, 1);
#ifdef DEBUG
		    cout << "XDet:\n" << XDet_itj[i][t][j] 
			 << "\ndetParmVec:\n" << detParmVec
			 << "\np:\n" << p << endl;
#endif
		    // get detection matrix for obs j
		    detMatObs = (*detMatPtr)(rbind(p, beta));

		    // pull out the proper column, according to y
		    detVecObs = detMatObs(_,y_itj[i][t][j]);
		} else {
		    // for missing data, obs prob mult by 1
		    detVecObs = ones(K + 1, 1);
		}
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
	    // allow user interrupts
	    R_CheckUserInterrupt();              
	}
    }

#ifdef DEBUG
    cout << "negll: " << negLogLike << endl;
#endif

    return negLogLike;
}

  
void HiddenMarkovModel::fitmodel(const Matrix<>& inits) {

    // set random number generator
    mersenne myrng;

    try{
	mle = scythe::BFGS(*this, inits, myrng, 1e5, 1e-5, false);
    } catch (scythe_convergence_error err) {
	cout << "BFGS did not converge." << endl;	    
    }

    hessian = scythe::hesscdif(*this, mle);

    const Matrix<> logitPhiParms = mle(nDP + K, 0, nP - 1, 0);
    const Matrix<> phiParms =  Hphi * plogis(logitPhiParms, 0, 1);
    phi = t((*phiMatPtr)(phiParms));
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
		  double *mle,
		  double *hessian,
		  double *negloglike,
		  double *phi)
    {

	Matrix<> Hdet = Matrix<> (*nDMP_un, *nDMP, Hdet_in);
	Matrix<> Hphi = Matrix<> (*nPhiP_un, *nPhiP, Hphi_in);
	Matrix<> inits(*nP, 1, true, 0);

	HiddenMarkovModel hmm(nDMP, nDCP, nDP, nDYP,
			      nSP, nPhiP, nP, nDMP_un, nPhiP_un,
			      K, M, J, nY, Hdet, Hphi, yearly_det);
	hmm.loadData(y_itj, XDet_itjk, ncXDet);

	if (*K == 2) {
	    cout << "entered 3 state" << endl;
	    phiMatrix3 phiMat;
	    detMatrix3 detMat;
	    phiMatrix *phiMatPtr = &phiMat;
	    detMatrix *detMatPtr = &detMat;
	    hmm.setMatrices(phiMatPtr, detMatPtr);
	} else {
	    cout << "model instantiated" << endl;
	    phiMatrix4 phiMat;
	    detMatrix4 detMat;
	    phiMatrix *phiMatPtr = &phiMat;
	    detMatrix *detMatPtr = &detMat;
	    hmm.setMatrices(phiMatPtr, detMatPtr);
	}

	// find the mle
	hmm.fitmodel(inits);

	// store the results
	for (int i = 0; i < *nP; ++i)
	    mle[i] = hmm.mle(i);
	for (int i = 0; i < (*nP)*(*nP); ++i)
	    hessian[i] = hmm.hessian(i);
	*negloglike = hmm(hmm.mle);
	for (int i = 0; i < (*K + 1) * (*K + 1); ++i)
	    phi[i] = hmm.phi(i);


    }
}
