///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_stat_TBeltLH__
#define __Stntuple_stat_TBeltLH__

#include "TNamed.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TRandom3.h"
#include "TGraph.h"

#include "Stntuple/stat/TBelt.hh"

namespace stntuple {
//-----------------------------------------------------------------------------
class TBeltLH: public TBelt {
public:

  struct Hist_t {
    TH1D*    fProb;
    TH1D*    fLlh;
    TH1D*    fLlhInterval;
    TH1D*    fInterval;
    TH1D*    fBeltLo;
    TH1D*    fBeltHi;
    TH1D*    fBeltNO;
    TH1D*    fBeltUL;
    TH1D*    fBeltNO1;
    // TH1D*    fSign[2];
    TGraph*  fCoverage;
  };

  Hist_t fHist;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
  TBeltLH(const char* Name, double CL = -1);

  virtual int   init_poisson_dist(double MuB, double MuS, int NObs = -1);

 // in presence of background, mu = mus+mub
  
  virtual int  construct_interval(double MuB, double MuS);
  virtual int  construct_belt    (double MuB, double SMin, double SMax, int NPoints, int NObs, double* P);

  virtual int  make_prob_hist();
  virtual void make_belt_hist();

  virtual int  test_coverage(double MuB, double SMin, double SMax, int NPoints);

  ClassDef(TBeltLH,0)
};

}
#endif
