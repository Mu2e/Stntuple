///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_stat_upper_limit__
#define __Stntuple_stat_upper_limit__

#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TGraph.h"

namespace stntuple {
//-----------------------------------------------------------------------------
class upper_limit: public TNamed {
public:

  enum {
    MaxNx =   50
  };

  struct Belt_t {
    double fSMin;
    double fSMax;
    double fDy;
    double fSign[MaxNx][2];      // fBelt.Sign[ix][0] = SMin[ix], fBelt.fSign[ix][1] = SMax[ix]
    
    int    fFillColor;
    int    fFillStyle;
    
    THStack* fStack;
    TH1D*    fLo;
    TH1D*    fHi;

    double   fXMin;
    double   fXMax;

    double   fYMin;
    double   fYMax;

    Belt_t();
    
    void  set_xmin(double X) { fXMin = X; }
    void  set_xmax(double X) { fXMax = X; }
    void  set_ymin(double Y) { fYMin = Y; }
    void  set_ymax(double Y) { fYMax = Y; }

    //    void  make_hist(int NObs);
    
  };

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
    TH1D*    fSign[2];
    TGraph*  fCoverage;
  };

  Hist_t fHist;
                                        // fBelt0 - unconditioned
  Belt_t fBelt0;
  Belt_t fBelt[MaxNx];

  int    fType;                      // default:0 if=1, try to force monotonic left edge
  double fCL;
  double fSump;

  double fMean;
  double fMuB;
  double fMuS;

  int    fNObs;
  
  int    fIxMin;
  int    fIxMax;

  long int fNExp;

  TRandom3 fRn;

  int    fIPMax;                    // index corresponding max(fProb);
  double fProb     [MaxNx];
  double fCProb    [MaxNx];
  double fFactorial[MaxNx];

  void SetNExp(long int NExp) { fNExp = NExp; }

  Belt_t* belt(int NObs) {
    if (NObs == -1) return &fBelt0;
    else            return &fBelt[NObs];
  }

  upper_limit(const char* Name, double CL = -1, int TYpe = 0);

                                        // 'N' - number of measured events
  
  int init_poisson_dist(double MuB, double MuS, int NObs = -1);
  
                                        // in presence of background, mu = mus+mub
  
  int construct_interval(double MuB, double MuS, int NObs = -1);
  int construct_belt    (double MuB, double SMin, double SMax, int NPoints,  int NObs = -1);

  int  make_prob_hist();
  
  void make_belt_hist(int NObs = -1);

  int test_coverage(double MuB, double SMin, double SMax, int NPoints);

  ClassDefOverride(upper_limit,0)
};

}
#endif
