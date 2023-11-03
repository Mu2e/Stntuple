///////////////////////////////////////////////////////////////////////////////
// fMode = 0: sample background momentum spectrum
//       = 1: sample signal     momentum spectrum
//
// initially, I forgot the fact that the signal and background
// means over the [103.6,104.9] are very close
// - the plan unchanged: normalize the weights properly,
//   assign each pseudoexperiment event the weight 
//   of the 2D probability density to get the values of pmean and sigmap
//   from the very beginning, consider the background as consisting of
//   two components: a uniform one and the DIO.
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_stat_TKinLH_hh__
#define __Stntuple_stat_TKinLH_hh__

#include "TRandom3.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TF1.h"
#include "THStack.h"

#include "TNamed.h"
#include "TFile.h"

#include "Stntuple/val/stntuple_val_functions.hh"
#include "Stntuple/stat/TBelt.hh"
//-----------------------------------------------------------------------------
namespace stntuple {
class TKinLH : public TBelt {
public:
                                        // Åxtracted histogram
  struct sdata {
    int    bin;
    double x;
  };
                                        // "experimental data"
  struct Data_t {
    int    fNEvents;
    double fP[100];                      // for the moment, assume 100 events is OK
  } fData;
  
  struct Debug_t {
    int      fRun;
    int      fConstructInterval;
    int      fConstructBelt;
    int      fTestCoverage;
    double   fMuMin;
    double   fMuMax;
  } fDebug;

  static int fDebug_QuickSort;

  struct Hist_t {
    TH1D*    fProb;
    TH1D*    fLlh;
    TH1D*    fLlhInterval;
    TH1D*    fInterval;
    THStack* fBelt;
    TH1D*    fBeltLo;
    TH1D*    fBeltHi;
    TH1D*    fBeltProb;
    TH1D*    fBeltLhdt;
    TGraph*  fCoverage;

    TH1F*    prob_sig;
    TH1F*    prob_bgr;
//-----------------------------------------------------------------------------
// fLoglhs[N][nb] is a histogram for N = Nb(bgr)+Ns(sig)
//-----------------------------------------------------------------------------
    TObjArray*      fLogLhs   [MaxNx];
    TObjArray*      fLogLhb   [MaxNx];
    TObjArray*      fLogLhr   [MaxNx];        // log(LHR) for given Nb and Ns (Nb+Ns=Ntot, Ntot fixed)
    TObjArray*      fLogLhrR  [MaxNx];        // log(LHR) / Ntot
    
    TH1D*           fLogLhrR_1[MaxNx];        // fLogLhrR for fiven ntot, summed over nb

    TH1D*           fLogLhrR_2[MaxNx];        // uniformly normalized, 2-sided distributions
                                              // sum of all for a given Ntot
    
    TH1D*           fSumLogLhrR_2;            // sum of all - for fiven MuB and MuS.

                                              // these are not used and are filled only in debug mode
    TH1D*    gen_pbgr;                        // generated momentum, background
    TH1D*    gen_psig;                        // generated momentum, signal
  };

  struct Interval_t {
    double fLlhrMin;
    double fLlhrMax;
    double fProbTot;
    double fPMax;
    int    fIMax;
  } fInterval;

  Hist_t   fHist;
  int      fColor;

  sdata*   fSortData;

  double   pmin;
  double   pmax;

  int      fInitialized;                // 1: histograms read in, 0: otherwise

  double   fMinLLHR;
  double   fMaxLLHR;
  
  TRandom3 fRng;
                                        // signal and background kinematic probability distributions
                                        // start with 1D momentum distributions, but, in general, those
                                        // are the 1D probability distributions integrated over all kinematic
                                        // variables.
  TF1*     fSig;
  TF1*     fBgr;

  void*    fTestHist;
                                        // momenta for testing calculations
  double   fTestMom1[10];
  double   fTestMom2[10];
  double   fTestMom3[10];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TKinLH(const char* Name, double CL = 0.9, double PMin=103.6, double PMax = 104.9, int Debug = 0);
  ~TKinLH();

  int    init();
                                        // signal parameterization with the integral \int f(x)*dx = 1
  
  double static f_sig(double* X, double * P);
  double static f_bgr(double* X, double * P);
  
  double lh_bgr (double P);
  double lh_sig (double P);
  double wt_data(double MuB, double MuS, int NObs, double* P, double* LlhrR);

  double bgr_mom();
  double sig_mom();

  int  construct_interval(double MuB, double MuS);

                                        // NObs is the number of observed events, P - array of their momenta
                                        // in principle, could pass a pointer (or a reference) to a vector
    
  virtual int  construct_belt    (double MuB, double SMin, double SMax, int NPoints, int NObs, double* P);
  virtual int  test_coverage     (double MuB, double SMin, double SMax, int NPoints);

  static int   partition(sdata arr[], int low , int high, double pivot);
  static void  quickSort(sdata arr[], int low , int high);

     // with the kinematic distributions in, NObs = -1 simply doesn't make sense
  
  int    run(int NObs, int NPe = 1000000);

  int    generate();

  int    save_hist(const char* Filename, const char* Opt = "append");
  int    read_hist(const char* Filename);
  
  void   plot_interval();
  void   make_belt_hist();

  virtual void Print(const char* Option_t) const;

  ClassDef(TKinLH,0)
};
}
#endif
