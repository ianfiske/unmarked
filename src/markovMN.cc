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
#include <R_ext/Utils.h> // includes R's interrupt capabilities
#include <R_ext/Print.h> // includes R's printing system
#include <vector>
#include <map>
#include <string>
#include "phimatrix.h"
#include "detmatrix.h"

//#define DEBUG2
//#define DEBUG

using namespace std;
using namespace scythe;
using std::vector;

static void InitializePhiMap();

// base class for the library
class HiddenMarkovModel {
public:
  double operator() (const Matrix<>& parms);
  unsigned int gradFit(const Matrix<> &detInits,
		       const Matrix<> &phiInits,
		       const Matrix<> &psiInits);
  //  unsigned int fitmodel(const Matrix<> &inits);
  void loadData(int *y_itj_in, double *XDet_itjk, const int *ncXDet);
  void setPhiMat(phiMatrix *phiMat);
  void setDetMat(detMatrix *detMat);
  //  inline Matrix<> getDetParms(const Matrix<> &parms);
  inline const Matrix<> getPhi(const Matrix<> &parms);
  inline const Matrix<> getPsi(Matrix<> psiParms);
  inline double forward( const Matrix<>& detParms,
			 const Matrix<>& phi, const Matrix<>& psi);
  inline void backward( const Matrix<>& detParms,
			const Matrix<>& phi, const Matrix<>& psi);
  void forwardBackward(const Matrix<>& psi, const Matrix<>& phi,
		       const Matrix<>& detParms);
  void getTransP(const Matrix<>& psi, const Matrix<>& phi,
		 const Matrix<>& detParms);
  Matrix<> getDetVec(const Matrix<>& detParmVec,
		     const int& i, const int& t) const;
  int EM(double tolerance, Matrix<> detInits,
	 Matrix<> phiInits, Matrix<> psiInits);
  Matrix<> getInits(const Matrix<> &starts);

  // constructor
  explicit HiddenMarkovModel(const int *nDMP,
			     const int *nDP,
			     const int *nSP,
			     const int *nPhiP,
			     const int *nP,
			     const int *nDP_un,
			     const int *nPhiP_un,
			     const int *K,
			     const int *M,
			     const int *J,
			     const int *nY,
			     const Matrix<> Hdet,
			     const Matrix<> Hphi,
			     const Matrix<> Hpsi,
			     const int *trace,
			     const int *homoDet)
    :Hdet(Hdet), Hphi(Hphi), Hpsi(Hpsi),
     nDMP(*nDMP), nDP(*nDP), nSP(*nSP),
     nPhiP(*nPhiP), nP(*nP), nDP_un(*nDP_un),
     nPhiP_un(*nPhiP_un),
     K(*K), M(*M), J(*J), nY(*nY), trace(*trace), homoDet(*homoDet)
  { }
  HiddenMarkovModel() {};

  Matrix<> Hdet;
  Matrix<> Hphi;
  Matrix<> Hpsi;
  int nDMP, nDP, nSP, nPhiP, nP, nDP_un, nPhiP_un,
    K, M, J;
  int nY;
  bool trace;
  bool homoDet;
  vector< vector< vector <unsigned int> > > y_itj;

  // Design matrix for each site, time, observation. dim:
  // K x ncol(XDet)
  vector< vector< vector < Matrix<> > > > XDet_itj;

  // forward-backward algorithm storage (notation of Rabiner 1989)
  vector< vector< Matrix<> > > forP;  // forward probs
  vector< vector< Matrix<> > > backP;   // backward probs
  vector< vector< Matrix<> > > smoothP;  // P(x_t = i | all y's and model)
  vector< vector< Matrix<> > > transP; // P(x_t = i, x_t+1 = j | O, lam)

  // Matrix<> mle; // will store the MLE after fitting and inits before
  Matrix<> hessian;


  // these store mle for EM algorithm
  Matrix<> detParms;
  Matrix<> phiParms;
  Matrix<> psiParms;
  Matrix<> phi;

  Matrix<> starts;

  phiMatrix *phiMatPtr;  // phi matrix functor pointer
  detMatrix *detMatPtr;  // det matrix functor pointer
  string phiMatString;
  friend class detLike;
  friend class phiLike;
};

const Matrix<> HiddenMarkovModel::getPhi(const Matrix<>& phiParms) {
  //  const Matrix<> logitPhiParms = parms(nDP + K, 0, nP - 1, 0);
  //Matrix<> phiParms, phi;
  /*  if (phiMatString.compare("cumlogit") != 0 ) {
    phiParms = Hphi * plogis(phiParms, 0, 1);
    phi = (*phiMatPtr)(phiParms);
    } else {*/

  return((*phiMatPtr)(Hphi * phiParms));
}

const Matrix<> HiddenMarkovModel::getPsi(Matrix<> psiParms) {
  Matrix<> Zero(0);
  psiParms = Hpsi * psiParms;
  psiParms = rbind(Zero, psiParms);
  Matrix<> psi = scythe::exp(psiParms);
  psi /= sum(psi);
  return(psi);
}

Matrix<> getBaselineLogit(const Matrix<>& psi) {
  Matrix<> baseLine(psi.rows() - 1, 1);
  for(int i = 0; i < psi.rows() - 1; i++) {
    baseLine[i] = log(psi[i]/psi[psi.rows()-1]);
  }
  return baseLine;
}


void HiddenMarkovModel::loadData(int *y_itj_in,
				 double *XDet_itjk,
				 const int *ncXDet)
{
  Matrix<> X;
  for (int i = 0; i < M; ++i) {
    y_itj.push_back(vector< vector< unsigned int > > () );
    XDet_itj.push_back(vector< vector< Matrix<> > > () );
    for (int t = 0; t < nY; ++t) {
      y_itj[i].push_back(vector< unsigned int > () );
      XDet_itj[i].push_back(vector< Matrix<> > () );
      for (int j = 0; j < J; ++j) {
	y_itj[i][t].push_back(*(y_itj_in++));

	// note that this is the tranpose of XDet
	X = Matrix<> (*ncXDet, nDMP, XDet_itjk);
	//

	XDet_itjk += (nDMP) * (*ncXDet);  // this line will change XDet changes.
	XDet_itj[i][t].push_back(scythe::t(X));
      }
    }
  }
}

double HiddenMarkovModel::forward(const Matrix<>& detParms,
				  const Matrix<>& phi,
				  const Matrix<>& psi) {

  Matrix<> detMatObs(K + 1, K + 1, false, 0);
  Matrix<> detVecObs;
  Matrix<> detVec;
  Matrix<> psiSite;
  double negLogLike = 0;

#ifdef DEBUG
  cout << "psi\n" << psi << endl;
  cout << "phiParms\n" << phiParms
       <<"\nphi\n" << phi << endl;
#endif

  if(homoDet) {
	  for(int k = 0; k < K + 1; ++k) {
		  detMatObs(_,k) = (*detMatPtr)(detParms, k);
	  }
  }

  for (int i = 0; i < M; ++i) {
      forP.push_back(vector< Matrix<> > () );
    // Initialize psi for site i.
    psiSite = psi;

    for (int t = 0; t < nY; ++t) {
      // initialize detVec to ones.
      // this will be the product the observation's detVec's.

      if(homoDet) { // if homogenous detection, pass constant det matrix instead of parms
    	  detVec = getDetVec(detMatObs, i, t);
      } else {
    	  detVec = getDetVec(detParms, i, t);
      }

      // mult by appropriate detMat and phi
      psiSite = psiSite % detVec;

      forP[i].push_back(psiSite);

      // if last time, then instead of multiplying by phi^t, just
      // sum up the vector.  this is the likelihood for site i.
      if (t < (nY - 1))
    	  psiSite = phi * psiSite;
      else
    	  negLogLike += -log(sum(psiSite));


      /*      if(isnan(negLogLike) | isinf(negLogLike)) {
	Rprintf("NLL is NAN.");
	abort();
	}*/


#ifdef DEBUG2
      cout << "detVec:\n" << detVec
	   << "psiSite:\n" << psiSite << endl;
#endif
      // allow user interrupts
      R_CheckUserInterrupt();
    }
#ifdef DEBUG2
    cout << "site  " << i << " nll: " << negLogLike << endl;
#endif
  }

  //    if (negLogLike == numeric_limits<double>::infinity()) negLogLike = 1e5;
#ifdef DEBUG
  //  cout << "negll: " << negLogLike << endl;
#endif
  return negLogLike;
}

void HiddenMarkovModel::backward(const Matrix<>& detParms,
				 const Matrix<>& phi, const Matrix<>& psi)

{
  Matrix<> detMatObs(K + 1, K + 1, false, 0);
  Matrix<> detVecObs;
  Matrix<> detVec;
  Matrix<> psiSite;
  Matrix<> backP;

  for (int i = 0; i < M ; ++i) {

    this->backP.push_back(vector< Matrix<> > () );

    backP = ones(K + 1, 1); // \beta_T = 1
    //	cout << "site "<< i <<"  pushed first back site vector" << endl;
    for (int t = nY - 1; t >= 0 ; t--) {

      // note that these betas are actually in reverse order
      this->backP[i].push_back(backP);

      detVec = getDetVec(detParms, i,t);

      backP = scythe::t(phi) * (detVec % backP);
    }
  }
}

Matrix<> HiddenMarkovModel::getDetVec(const Matrix<>& detParms,
				      const int& i,
				      const int& t) const
{
	Matrix<> mp, detMatObs, detVecObs;
	Matrix<> detVec = ones(K + 1, 1);
	for (int j = 0; j < J; ++j) {
		if (y_itj[i][t][j] != 99) {

			if(homoDet) {  // note the abuse of notation... detParms is actually detMatObs in this case. (ugly code!)
				detVec = detVec % detParms(_,y_itj[i][t][j]);
			} else {
				// compute matrix parameters
				mp = XDet_itj[i][t][j] * detParms;

				detVecObs = (*detMatPtr)(mp, y_itj[i][t][j]);


				// update yearly detection vector
				detVec = detVec % detVecObs;
			}
		}
	}

  return detVec;
}

void HiddenMarkovModel::setPhiMat(phiMatrix *phiMatPtr)
{
  this->phiMatPtr = phiMatPtr;
}

void HiddenMarkovModel::setDetMat(detMatrix *detMatPtr)
{
  this->detMatPtr = detMatPtr;
}


double HiddenMarkovModel::operator() (const Matrix<>& parms)
{
  forP.clear();
  backP.clear();
  transP.clear();
  smoothP.clear();

  const Matrix<> psi = getPsi(parms(0, 0, nSP - 1, 0));
  const Matrix<> phi = (*phiMatPtr) (Hphi * parms(nSP, 0, nSP - 1 + nPhiP, 0));
  const Matrix<> detParms = Hdet * parms(nSP + nPhiP, 0,
					 nSP - 1 + nDP + nPhiP, 0);

  double negLogLike = forward(detParms, phi, psi);
  R_CheckUserInterrupt();
  return(negLogLike);
}


void HiddenMarkovModel::forwardBackward(const Matrix<>& psi,
					const Matrix<>& phi,
					const Matrix<>& detParms) {

  forward(detParms, phi, psi);
  backward(detParms, phi, psi);
  Matrix<> gamma, fP, bP;

  for(int i = 0; i < M; i++) {
    smoothP.push_back(vector< Matrix<> > () );
    for (int t = 0; t < nY; t++) {
      fP = (this->forP)[i][t];
      bP = (this->backP)[i][nY - 1 - t]; //these were in reverse order
      gamma = 1/(scythe::t(fP) * bP) * (fP % bP);
      smoothP[i].push_back(gamma);
      /*      for(int k=0; k <= K; k++) {
	if (isnan(gamma(k)))
	  cout << "i " << i << " t " << t << "\n\n" << gamma << endl;;
	  }*/
    }
  }
  //    cout << "forward-backward complete" << endl;
}

// computes the probability of each transition
void HiddenMarkovModel::getTransP(const Matrix<>& psi,
				  const Matrix<>& phi,
				  const Matrix<>& detParms)
{

  Matrix<> tP(K + 1, K + 1, false);
  Matrix<> detVec;

  for(int i = 0; i < M; i++) {

    transP.push_back(vector< Matrix<> > () );
    for(int t = 0; t < nY - 1; t++) {

      detVec = getDetVec(detParms, i, t + 1);

      for(unsigned int k = 0; k <= K; k++) {
	for(unsigned int l = 0; l <= K; l++) {
	  // recall that phi is still transposed (col-stoc) => tP is col-stoc
	  // recall bP is backwards
	  tP(l,k) = forP[i][t](k) * phi(l,k) * detVec(l) * backP[i][nY-t-2](l);
	}
      }

      transP[i].push_back(tP/sum(tP));
      //      cout << "tp\n" << transP[i][t] << "sp\n" << smoothP[i][t] << endl;
    }
  }


}

class detLike {
public:
  double operator()(const Matrix<>& detParms_con) {

    const Matrix<> detParms_un = hmm.Hdet * detParms_con;
    Matrix<> detVec;
    Matrix<> sP(hmm.K + 1, 1, false);

    double negloglike = 0;
    for (int i = 0; i < hmm.M; i++) {
      for (int t = 0; t < hmm.nY; t++) {
	detVec = hmm.getDetVec(detParms_un, i, t);
	sP = hmm.smoothP[i][t];
	for (int k = 0; k <= hmm.K; k++) {
	  //detvec_k should only be zero when sp_k is also 0
	  //so the "if" takes care of inf*0 = 0
	  if(detVec[k] != 0)
	    negloglike -= log(detVec[k]) * sP[k];
	}
      }
    }
    R_CheckUserInterrupt();
    return(negloglike);
  }

  detLike(HiddenMarkovModel H): hmm(H) { };

  HiddenMarkovModel hmm;
};

class phiLike {
public:
  double operator() (const Matrix<>& phiParms_con) {
    const Matrix<> phiParms = hmm.Hphi * phiParms_con;
    const Matrix<> phi = (*hmm.phiMatPtr)(phiParms);

    double phinll = 0;
    for(int i = 0; i < hmm.M; i++) {
      for(int t = 0; t < hmm.nY - 1; t++) {
	for(int k = 0; k <= hmm.K; k++) {
	  for(int l = 0; l <= hmm.K; l++) {
	    phinll -= hmm.transP[i][t](l,k) * log(phi(l,k));
	  }
	}
      }
    }
    R_CheckUserInterrupt();
    return(phinll);

  }

  phiLike(HiddenMarkovModel H): hmm(H) { };
  HiddenMarkovModel hmm;

};

// EM IS BROKEN WITH CONSTRAINED PSI!  RELIMPMLEMENT M-STEP FOR PSI!
int HiddenMarkovModel::EM(double tolerance, Matrix<> detInits,
			  Matrix<> phiInits, Matrix<> psiInits) {

  Matrix<> detParms_con_next = detInits;
  Matrix<> phiParms_con_next = phiInits;
  Matrix<> psi_next = getPsi(psiInits);

  double iterDiff;
  mersenne myrng;
  detLike dlike(*this);
  phiLike plike(*this);

  Matrix<> phi, psi, detParms, denom, detParms_con,
    detParms_un, detParms_next,
    phiParms_con;

  Matrix<> zeroVec(K + 1, 1, true, 0);
  Matrix<> zeroMat(K + 1, K + 1, true, 0);

  do {
    iterDiff = 0;

    psi = psi_next;
    phiParms_con = phiParms_con_next;
    detParms_con = detParms_con_next;

    // E-step - smooth to get forP, backP, smoothP, and transP
    phi = (*phiMatPtr)(Hphi * phiParms_con);
    forP.clear();
    backP.clear();
    transP.clear();
    smoothP.clear();
    detParms_un = Hdet * detParms_con;
    forwardBackward(psi, phi, detParms_un);
    getTransP(psi, phi, detParms_un);

    // get likelihood for debugging
    double nll = 0;
    for (int i = 0; i < M; i++) {
      nll -= log(sum(forP[i][nY - 1]));
    }
    cout << "nll = " << nll << endl;


    // M-step
    // update psi;
    psi_next = zeroVec;
    for (int i = 0; i < M; i++) {
      psi_next += smoothP[i][0];
    }
    psi_next = psi_next / M;
    iterDiff += sum((psi - psi_next)^2);

    // update phi;
    plike.hmm.transP = transP;
    phiParms_con_next = scythe::BFGS(plike, phiParms_con,// + phiParmDev,
				     myrng, 200, 1e-5, trace);
    iterDiff += sum((phiParms_con - phiParms_con_next)^2);

    // update detparms numerically;
    dlike.hmm.smoothP = smoothP;
    detParms_con_next = scythe::BFGS(dlike, detParms_con,// + detParmDev,
				     myrng, 200, 1e-5, trace);
    iterDiff += sum((detParms_con_next - detParms_con)^2);

    cout << "iterDiff = " << iterDiff << endl;
    R_CheckUserInterrupt();
  } while (iterDiff > tolerance);

  this->phiParms = phiParms_con_next;
  this->detParms = detParms_con_next;
  //  this->psi = psi_next;
  this->phi = scythe::t((*phiMatPtr)(Hphi * phiParms_con_next));
  // Matrix<> mle = getBaselineLogit(this->psi);
  // mle = rbind(mle, this->phiParms);
  // mle = rbind(mle, this->detParms);
  //  this->hessian = scythe::hesscdif(*this, mle);

  return 1;
}


unsigned int HiddenMarkovModel::gradFit(const Matrix<> &detInits,
					const Matrix<> &phiInits,
					const Matrix<> &psiInits)
{
  mersenne myrng;
  Matrix<> inits = rbind(psiInits, phiInits);
  inits = rbind(inits, detInits);
  Matrix<> mle(nP, 1);
  try
    {
      Rprintf("Entering BFGS...\n");
      mle = scythe::BFGS(*this, inits, myrng, 200, 1e-8, trace);
    }
  catch (scythe_convergence_error err)
    {
      Rprintf("BFGS did not converge.\n");
      return(0);
    }

//TODO why does convergence failure sometimes appear to not get caught?

 /*
 * The scythe_convergence_error is not caught, but yet nonsensical fits get returned
 * with a hessian full of NaN's.
 *
 * Maybe add manual checking of hessian?  do this in R?
 *
 */

  Rprintf("BFGS converged.\n");
  hessian = scythe::hesscdif(*this, mle);
  psiParms = mle(0, 0, nSP - 1, 0);
  phiParms =  mle(nSP, 0, nSP - 1 + nPhiP, 0);
  detParms = mle(nSP + nPhiP, 0, nSP - 1 + nDP + nPhiP, 0);
  phi = scythe::t((*phiMatPtr)(Hphi * phiParms));
  Matrix<> psi = getPsi(psiParms);
  forwardBackward(psi, phi, Hdet * detParms);

  return(1);
}

Matrix<> HiddenMarkovModel::getInits(const Matrix<>& starts) {
  Rprintf("Getting initial values...\n");
  double delta = 1, nll = (*this)(starts), nll_try;
  Matrix<> parms = starts;
  Matrix<> parms_try;
  while(delta > 1e-2) {
    for(int i = 0; i < nP; i++) {
      parms_try = parms;
      parms_try[i] -= delta;
      nll_try = (*this)(parms_try);
      if (nll_try < nll) {
	nll = nll_try;
	parms = parms_try;
      }

      parms_try[i] += 2*delta;
      nll_try = (*this)(parms_try);
      if (nll_try < nll) {
	nll = nll_try;
	parms = parms_try;
      }
    }
    delta /= 1.2;
  }

  return(parms);
}

extern "C" {
  /*  void phi4state(const double* phiParms,
		 double* phi) {
    phiMatrix4 phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (12, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 16; ++i)
      phi[i] = p(i);
      }

  void phi4state4(const double* phiParms,
		  double* phi) {
    phiMatrix4p4 phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (4, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 16; ++i)
      phi[i] = p(i);
      }

  void phi3state(const double* phiParms,
		 double* phi) {
    phiMatrix3 phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (6, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 9; ++i)
      phi[i] = p(i);
  }

  void phi2state(const double* phiParms,
		 double* phi) {
    phiMatrix2 phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (2, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 4; ++i)
      phi[i] = p(i);
  }


*/

  void phi4stateAR(const double* phiParms,
		   double* phi) {
    phiMatrix4ar phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (10, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 16; ++i)
      phi[i] = p(i);
  }


  void phi4logit(const double* phiParms,
		 double* phi) {

    phiLogit4 phiMatFun;
    Matrix<> phiParmsMat = Matrix<> (12, 1, phiParms);
    Matrix<> p = phiMatFun(phiParmsMat);

    for (int i = 0; i < 4; ++i)
      phi[i] = p(i);
  }


  void findMLE(int *y_itj,
	       double *XDet_itjk,
	       const int *ncXDet,
	       const int *nDMP,
	       const int *nDP,
	       const int *nSP,
	       const int *nSP_un,
	       const int *nPhiP,
	       const int *nP,
	       const int *nDP_un,
	       const int *nPhiP_un,
	       const double *Hdet_in,
	       const double *Hphi_in,
	       const double *Hpsi_in,
	       const int *K,
	       const int *M,
	       const int *J,
	       const int *nY,
	       double *negloglike,
	       double *phi,
	       double *psiParms,
	       double *detParms,
	       double *phiParms,
	       double *smooth,
	       const int *smoothOnly,
	       const char** phiTypeCString,
	       double *hessian,
	       const int *EM,
	       double *psiInit,
	       double *phiParmInit,
	       double *detParmInit,
	       int *seekInits,
	       int *trace,
	       int *homoDet,
	       int *arDet,
	       int *convergence)
  {

    Matrix<> Hdet = Matrix<> (*nDP_un, *nDP, Hdet_in);
    Matrix<> Hphi = Matrix<> (*nPhiP_un, *nPhiP, Hphi_in);
    Matrix<> Hpsi = Matrix<> (*nSP_un, *nSP, Hpsi_in);
    Matrix<> detInits = Matrix<> (*nDP, 1, detParmInit);
    Matrix<> phiInits = Matrix<> (*nPhiP, 1, phiParmInit);
    Matrix<> psiInits = Matrix<> (*nSP, 1, psiInit);


    HiddenMarkovModel hmm(nDMP, nDP,
			  nSP, nPhiP, nP, nDP_un, nPhiP_un,
			  K, M, J, nY, Hdet, Hphi, Hpsi, trace, homoDet);
    hmm.loadData(y_itj, XDet_itjk, ncXDet);

    string phiTypeString(*phiTypeCString);
    hmm.phiMatString = phiTypeString;
    InitializePhiMap();


    // phiMatrix2 phiMat2;
    // phiMatrix *phiMatPtr2 = &phiMat2;
    // phiMatrix3 phiMat3;
    // phiMatrix *phiMatPtr3 = &phiMat3;
    // phiMatrix4 phiMat4;
    // phiMatrix *phiMatPtr4 = &phiMat4;
    phiMatrix4ar phiMat4ar;
    phiMatrix *phiMatPtr4ar = &phiMat4ar;
    //    phiMatrix4p4 phiMat4p4;
    //    phiMatrix *phiMatPtr4p4 = &phiMat4p4;
    //    cumLogitPhiMatrix phiMatCL;
    //    phiMatCL.K = *K;
    //    phiMatrix *phiMatPtrCL = &phiMatCL;
    phiLogit2 phiMatL2;
    phiMatrix *phiMatPtrL2 = &phiMatL2;
    phiLogit3 phiMatL3;
    phiMatrix *phiMatPtrL3 = &phiMatL3;
    phiLogit4 phiMatL4;
    phiMatrix *phiMatPtrL4 = &phiMatL4;
    phiLogit4ar phiMatL4ar;
    phiMatrix *phiMatPtrL4ar = &phiMatL4ar;
    phi4random phiMat4R;
    phiMatrix *phiMatPtr4R = &phiMat4R;


    // set the Phi Matrix
    switch (phi_mapStringValues[phiTypeString]) {
    // case evMacKenzie:
    //   {
    // 	hmm.setPhiMat(phiMatPtr2);
    // 	if(*smoothOnly == 0) Rprintf("Fitting MacKenzie model ...");
    //   }
    //   break;
    // case evThreeState:
    //   {
    // 	hmm.setPhiMat(phiMatPtr3);
    //   }
    //   break;
    // case evFourState:
    //   {
    // 	hmm.setPhiMat(phiMatPtr4);
    // 	if(*smoothOnly == 0) Rprintf("Fitting 4 state model ...");
    //   }
    //   break;
    case evFourStateAR:
      {
	hmm.setPhiMat(phiMatPtr4ar);
        Rprintf("Fitting 4 state AR model ...");
      }
      break;
    case evLogit4:
      {
	hmm.setPhiMat(phiMatPtrL4);
	Rprintf("Fitting 4 state logit model ... \n");
      }
      break;
    case evLogit3:
      {
	hmm.setPhiMat(phiMatPtrL3);
	Rprintf("Fitting 3 state logit model ...\n");
      }
      break;
    case evLogit2:
      {
	hmm.setPhiMat(phiMatPtrL2);
	Rprintf("Fitting 2 state logit model ...\n");
      }
      break;
    case evLogit4ar:
      {
	hmm.setPhiMat(phiMatPtrL4ar);
      }
      break;
    case evFourStateRandom:
      {
	hmm.setPhiMat(phiMatPtr4R);
      }
      break;
    default:
      Rprintf("Error: unknown phi matrix specified.");
      return;
      break;
    }

    //detMatrix2 detMat2;
    //detMatrix *detMatPtr2 = &detMat2;
    detLogit3 detMat3;
    detMatrix *detMatPtr3 = &detMat3;
    detLogit2 detMat2;
    detMatrix *detMatPtr2 = &detMat2;
    detLogit4 detMat4;
    detMatrix *detMatPtr4 = &detMat4;
    detLogit4ar detMat4ar;
    detMatrix *detMatPtr4ar = &detMat4ar;
    // set the detection matrix
    switch (*K) {
    case 1:
      {
	hmm.setDetMat(detMatPtr2);
      }
      break;
    case 2:
      {
    	hmm.setDetMat(detMatPtr3);
      }
      break;
    case 3:
      {
	if(*arDet == 1)
	  hmm.setDetMat(detMatPtr4ar);
	else
	  hmm.setDetMat(detMatPtr4);
      }
      break;
    default:
      Rprintf("There is no detection matrix for K = %u.",*K);
      break;
    }

    unsigned int fitSuccess;

    if(*seekInits == 1){
      Matrix<> starts = rbind(psiInits, phiInits);
      starts = rbind(starts, detInits);
      Matrix<> inits = hmm.getInits(starts);
      psiInits = inits(0, 0, *nSP - 1, 0);
      phiInits = inits(*nSP, 0, *nSP - 1 + *nPhiP, 0);
      detInits = inits(*nSP + *nPhiP, 0, *nSP - 1 + *nDP + *nPhiP, 0);
    }

    if(*EM == 1) {
      fitSuccess = hmm.EM(1e-3, detInits, phiInits, psiInits);
    } else {
      fitSuccess = hmm.gradFit(detInits, phiInits, psiInits);
    }


    if (fitSuccess == 0) {
      Rprintf("Error. Fit did not converge.\n");
      *negloglike = numeric_limits<double>::infinity();
      *convergence = 1;
      return;
    }

    *convergence = 0;
    Rprintf("Storing results...\n");

    *negloglike = 0;
    for (int i = 0; i < *M; i++) {
      *negloglike -= log(sum(hmm.forP[i][*nY - 1]));
    }

    for (int i = 0; i < *nDP; i++)
      detParms[i] = hmm.detParms(i);
    for (int i = 0; i < *nPhiP; i++)
      phiParms[i] = hmm.phiParms(i);
    for (int i = 0; i < (*K + 1) * (*K + 1); ++i)
      phi[i] = hmm.phi(i);
    for (int i = 0; i <= *nSP; ++i)
      psiParms[i] = hmm.psiParms(i);
    for (int i = 0; i < (*nP)*(*nP); ++i)
      hessian[i] = hmm.hessian(i);

    unsigned int q = 0;
    for (int i = 0; i < *M; i++) {
      for (int t = 0; t < *nY; t++) {
	for (int k = 0; k <= *K; k++, q++) {
	  smooth[q] = hmm.smoothP[i][t](k);
	  //		    cout << q << " " << smooth[q] << endl;
	}
      }
    }

    //	cout << "end." << endl;
  }


}

void InitializePhiMap()
{
  phi_mapStringValues["2state"] = evMacKenzie;
  phi_mapStringValues["3state"] = evThreeState;
  phi_mapStringValues["4state"] = evFourState;
  phi_mapStringValues["4state4"] = evFourState4;
  phi_mapStringValues["4stateAR"] = evFourStateAR;
  phi_mapStringValues["4random"] = evFourStateRandom;
  phi_mapStringValues["cumlogit"] = evCumLogit;
  phi_mapStringValues["logit4"] = evLogit4;
  phi_mapStringValues["logit4ar"] = evLogit4ar;
  phi_mapStringValues["logit3"] = evLogit3;
  phi_mapStringValues["logit2"] = evLogit2;
  phi_mapStringValues["end"] = evEnd;
}
