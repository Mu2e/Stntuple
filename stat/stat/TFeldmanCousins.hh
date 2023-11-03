//
// if needed, modify interface, to accomodate for the model with systematics
//
#ifndef __Stntuple_stat_TFeldmanCousins__
#define __Stntuple_stat_TFeldmanCousins__

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TMath.h"
#include "TNamed.h"
#include  "TGraph.h"

#include "Stntuple/stat/model_t.hh"

namespace stntuple {
class TFeldmanCousins : public TNamed {
public:
  enum {
    MaxNx =   50,                   // max Poisson bin
    MaxNy = 1001		    // max steps in Mu
  };   

  struct Hist_t {
    TH1D*    fBgProb;
    TH1D*    fBsProb;
    TH1D*    fProb;                     // the same as fBsProb, but shifted by half-bin
    THStack* fBelt;
    TH1D*    fBeltLo;
    TH1D*    fBeltHi;
    TGraph*  fCoverage;
  } fHist;

  struct Belt_t {
    int    fNy;                   // use not to reinitialize
    int    fNObs;
    double fBgr;
    double fSMin;
    double fSMax;
    double fDy;
    double fSign[MaxNx][2];      // fSBelt[ix][0] = SMin[ix], fSBelt[ix][1] = SMax[ix]
    int    fFillColor;           // for belt histogram
    int    fFillStyle;
  } fBelt;

  struct DebugLevel_t {
    int    fAll;
    int    fConstructBelt;
    int    fConstructInterval;
    int    fUpperLimit;
    int    fTestCoverage;
    double fMuMin;
    double fMuMax;
  } fDebug;
  
  double   fCL;
  double   fLog1mCL;			// log(1-fCL)
  double   fMuB;
  double   fMuS;
  int      fNObs;
  
  TRandom3 fRn;

  double   fBestSig  [MaxNx];
  double   fBestProb [MaxNx];

  double   fLhRatio  [MaxNx];
  int      fRank     [MaxNx];

  double   fBgProb   [MaxNx];		//
  double   fBsProb   [MaxNx];

  // double   fCumBgProb[MaxNx];           // fCumBgProb[i]: P(X<=i) for BGR-only hyp
  // double   fCumBsProb[MaxNx];

  double   fFactorial[MaxNx];	        // precalculate to speed things up
  
  // X(Bg) - Poisson random number sampled from Bg     distribution
  // X(Bs) - Poisson random number sampled from Bg+Sig distribution
  
  long int fNExp;	      // N(pseudo experiments) to throw

  int      fIxMin;
  int      fIxMax;
  int      fNSummed;          // likely, fIMax-fIMin+1
  double   fProb;

  int      fType;             // by default: 0 (1: biased observations)
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TFeldmanCousins(const char* Name,
		   double      CL,
		   int         DebugLevel = 0);
  
  ~TFeldmanCousins();
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void   SetCL          (double CL);

  void   SetNExp        (long int N) { fNExp =    N; };
  void   SetType        (int   Type) { fType = Type; };
  
  void   SetDebugLevel  (int Level ) { fDebug.fAll = Level; };
  
  int    ConstructInterval(double Bgr, double Sig, int NObs = -1);

  int    ConstructInterval(model_t* Model);
  
  int    ConstructBelt    (double Bgr, double SMin, double SMax, int NPoints, int NObs = -1);

  double Factorial(int Ix) { return fFactorial[Ix]; }

  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to prob=50%

  void   DiscoveryProb    (double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob);
  
  // plot discovery probability for a given background and signal range
  // calling makes sense only if CL=-1
  // discovery corresponds to NSig = 5
  // output: arrays of MuS and NSig, NPoints in size each
  // if NPoints=1, calculate NSig only for SMin
  
  void   DiscoveryProbMean(double   MuB  , double SMin, double SMax, int NPoints, double* MuS, double* NSig);
  void   DiscoveryProbMean(model_t* Model, double SMin, double SMax, int NPoints, double* MuS, double* NSig);

  // MuS : signal corresponding to 50% prob
  // Prob: probability corresponding to the signal
  // (should be close to 50%)
  
  void   DiscoveryMedian  (double   MuB  , double SMin, double SMax, double* MuS, double* Prob) ;
  void   DiscoveryMedian  (model_t* Model, double SMin, double SMax, double* MuS, double* Prob) ;

  // void   Init           (double Bgr, double Sig);

  void   InitPoissonDist(double MuB, double MuS, double* Prob, int NObs = -1);

  void   MakeBeltHist();
  void   MakeProbHist();

  void   PrintData(const char* Title, char DataType, void* Data, int MaxInd);
  void   PrintProbs(int N);

  void   PlotDiscoveryProbMean(double MuB, double Mu2);

  int    SolveFor(double Val, const double* X, const double* Y, int NPoints, double* XVal);

  int    TestCoverage(double MuB, double SMin, double SMax, int NPoints, int NObs = -1);
    
  void   UpperLimit(double   MuB  , double SMin, double SMax, int NPoints, double* S, double* Prob);

  int    UpperLimit(double   MuB  , double SMin, double SMax, double* S, double* Prob);
  int    UpperLimit(model_t* Model, double SMin, double SMax, double* S, double* Prob);

  ClassDef(stntuple::TFeldmanCousins,0)
};
};
#endif
