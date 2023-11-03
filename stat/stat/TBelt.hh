///////////////////////////////////////////////////////////////////////////////
// base class for confidence interval / confidence belt studies
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_stat_TBelt__
#define __Stntuple_stat_TBelt__

#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TGraph.h"

namespace stntuple {
//-----------------------------------------------------------------------------
class TBelt: public TNamed {
public:

  enum {
    MaxNx =   50
  };

  struct Belt_t {
    double  fSMin;
    double  fSMax;
    double  fDy;
    double  fSign[MaxNx][2];      // fBelt.Sign[ix][0] = SMin[ix], fBelt.fSign[ix][1] = SMax[ix]
    
    int     fFillColor;
    int     fFillStyle;
    double  fXMin;
    double  fXMax;
    double  fYMin;
    double  fYMax;
    double* fLlhInterval;
    int     fLlhNPoints;
  } fBelt;

  struct Debug_t {
    int    fConstructBelt;
    int    fUpperLimit;
    int    fTestCoverage;
    double fMuMin;
    double fMuMax;
  } fDebug;
  
  struct Hist_t {
    TH1D*    fProb;
    TH1D*    fLh;
    TH1D*    fInterval;
    THStack* fBelt;
    TH1D*    fBeltLo;
    TH1D*    fBeltZr;
    TH1D*    fBeltHi;
    TH1D*    fSign[2];
    TGraph*  fCoverage;
  };

  Hist_t fHist;

  double fCL;
  double fSump;

  double fMean;
  double fMuB;
  double fMuS;

  int    fNObs;                         // N(observed events)
  
  int    fIxMin;
  int    fIxMax;

  long int fNExp;

  TRandom3 fRn;

  int    fIPMax;                    // index corresponding max(fProb);
  int    fRank     [MaxNx];
  double fProb     [MaxNx];
  double fCProb    [MaxNx];
  double fFactorial[MaxNx];

  void SetNExp(long int NExp) { fNExp = NExp; }

  void set_belt_xmax(double XMax) { fBelt.fXMax = XMax; }
  void set_belt_ymax(double YMax) { fBelt.fYMax = YMax; }

  TBelt(const char* Name, double CL = -1);

                                        // 'N' - number of measured events - both 'truncated' functions need to be removed....
                                        // 1. Zech
  int          init_truncated_poisson_dist(double MuB, int   NObs, double* Prob    );

                                        // 2. Mur
  int          init_truncated_poisson_dist(double MuB, double MuS, int NObs, double* Prob);

  int          init_poisson_dist          (double MuB, double* Prob);
  
  virtual int  init_poisson_dist          (double MuB, double MuS, int    NObs = -1);
  
                                        // in presence of background, mu = mus+mub
  
  virtual int  construct_interval(double MuB, double MuS);
  virtual int  construct_belt    (double MuB, double SMin, double SMax, int NPoints, int NObs, double* P);

  virtual int  make_prob_hist();
  virtual void make_belt_hist();

  virtual int  test_coverage(double MuB, double SMin, double SMax, int NPoints);

  ClassDef(TBelt,0)
};

}
#endif
