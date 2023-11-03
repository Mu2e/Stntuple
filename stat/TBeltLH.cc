///////////////////////////////////////////////////////////////////////////////
//
// fDebuglevel.fTestCoverage = 1: print all
// fDebuglevel.fTestCoverage = 2: print missed
///////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "Stntuple/stat/TBeltLH.hh"

ClassImp(stntuple::TBeltLH)

namespace stntuple {
//-----------------------------------------------------------------------------
TBeltLH::TBeltLH(const char* Name, double CL) : TBelt(Name,CL) {
  fHist.fProb       = nullptr;
  fHist.fLlh        = nullptr;
  fHist.fInterval   = nullptr;
  fHist.fBeltLo     = nullptr;
  fHist.fBeltHi     = nullptr;
  fHist.fBeltNO     = nullptr;
  fHist.fBeltUL     = nullptr;
  fHist.fCoverage   = nullptr;
}

//-----------------------------------------------------------------------------
int TBeltLH::init_poisson_dist(double MuB, double MuS, int NObs) {
  int rc = TBelt::init_poisson_dist(MuB, MuS, NObs);
  // in addition to initializing the Poisson distribution, rank its bins assuming the distribution
  // is falling on both sides from the maximum
  // ranking is unique for both sides
  // fIPMax - bin in N-space where the probability density reaches its maximum

  fRank[0] = fIPMax;

  int nranked = 1;
  
  int imax = fIPMax+1;
  int imin = fIPMax-1;

  while (nranked < MaxNx) {
    if ((imin >= 0) and (imax <=MaxNx)) {
      if (fProb[imax] > fProb[imin]) {
        fRank[nranked] = imax++;
      }
      else {
        fRank[nranked] = imin--;
      }
    }
    else {
      if (imax <= MaxNx) {
        fRank[nranked] = imax++;
      }
      else {
        fRank[nranked] = imin--;
      }
    }
    nranked++;
  }

  return rc;
}
//-----------------------------------------------------------------------------
// three parameters to maintain uniform interface
//-----------------------------------------------------------------------------
int TBeltLH::construct_interval(double MuB, double MuS) {
  // defined fIxMin, fIxMax, and fSumP
  
  int rc = init_poisson_dist(MuB,MuS);
  if (rc < 0) return rc;

  rc       = -1;
  fSump    = 0;
  fIxMin   = fIPMax;
  fIxMax   = fIPMax;

  for (int ix=0; ix<MaxNx; ix++) {
    int loc = fRank[ix];
    
    if      (loc < fIxMin) fIxMin = loc;
    else if (loc > fIxMax) fIxMax = loc;
    
    fSump += fProb[loc];
    if (fSump >= fCL) {
      rc       = 0;
      break;
    }
  }

  if (rc < 0) return rc;
//-----------------------------------------------------------------------------
// Crow-Gardner: see, if can move to the right to minimize the overcoverage
// fSump : covered by the interval
// fIxMin and IxMax: two ends of the interval in the N-space
//-----------------------------------------------------------------------------
  while (fIxMax < MaxNx) {
    double sp = fSump+ fProb[fIxMax+1]-fProb[fIxMin];
    if (sp > fCL) {
      fSump  = sp;
      fIxMin = fIxMin+1;
      fIxMax = fIxMax+1;
    }
    else {
      break;
    }
  }
  return rc;
}

//-----------------------------------------------------------------------------
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
//-----------------------------------------------------------------------------
  int TBeltLH::construct_belt(double MuB, double SMin, double SMax, int NPoints, int NObs, double* P) {

  fMuB = MuB;
  
  fBelt.fSMin = SMin;
  fBelt.fSMax = SMax;
  fBelt.fDy   = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 1;

  if ((fBelt.fLlhInterval) and (fBelt.fLlhNPoints = NPoints)) delete fBelt.fLlhInterval;

  fBelt.fLlhNPoints  = NPoints;
  fBelt.fLlhInterval = new double[5*NPoints];

  for (int i=0; i<5*NPoints; i++) {
    fBelt.fLlhInterval[i] = 0;
  }

 
  for (int i=0; i<NPoints; i++) {
    double mus = SMin+i*fBelt.fDy;

    int rc       = construct_interval(MuB,mus);
    double lhmax = fProb[fIPMax];

    if (rc == 0) {
      double llh_lo = -log(fProb[fIxMin])+log(lhmax);
      if (fIxMin < fIPMax) llh_lo = -llh_lo;

      double llh_hi = -log(fProb[fIxMax])+log(lhmax);
      
      fBelt.fLlhInterval[5*i  ] = llh_lo;
      fBelt.fLlhInterval[5*i+1] = llh_hi;

      double llh_nobs = 0;
      if (NObs >= 0) {
        llh_nobs = -log(fProb[NObs])+log(lhmax);
        if (NObs < fIPMax) llh_nobs = -llh_nobs;
      }
      fBelt.fLlhInterval[5*i+2] = llh_nobs;

      if (mus+MuB < NObs) fBelt.fLlhInterval[5*i+3] = -log(fProb[fIxMax])+log(lhmax);
      else                fBelt.fLlhInterval[5*i+3] = -log(fProb[fIxMin])+log(lhmax);

      fBelt.fLlhInterval[5*i+4] = -log(fProb[NObs])+log(lhmax);
    }
  }

  if (fDebug.fConstructBelt > 0) {
    for (int ix=0; ix<MaxNx; ix++) {
      printf("%32s %3i %12.5f %12.5f\n","",ix,fBelt.fSign[ix][0],fBelt.fSign[ix][1]);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TBeltLH::make_prob_hist() {

  if (fHist.fProb    ) delete fHist.fProb;
  if (fHist.fInterval) delete fHist.fInterval;
  
  TString title = Form("prob vs N, muB:%5.2f muS:%7.4f CL:%5.3f Nobs:%3i",fMuB,fMuS,fCL,fNObs);
  fHist.fProb     = new TH1D(Form("h_prob_%s",GetName()),title.Data(),MaxNx,-0.5,MaxNx-0.5);

  title = Form("Interval muB:%5.2f muS:%7.4f CL:%5.3f Nobs:%3i",fMuB,fMuS,fCL,fNObs);
  fHist.fInterval = new TH1D(Form("h_prob_i_%s",GetName()),title.Data(),MaxNx,-0.5,MaxNx-0.5);

  fHist.fInterval->SetFillStyle(3004);
  fHist.fInterval->SetFillColor(kRed+2);
  fHist.fInterval->SetLineColor(kRed+2);

  fHist.fLlh         = new TH1D(Form("h_llh_%s"  ,GetName()),"log(p/pmax)"         ,10000,-50,50);
  fHist.fLlhInterval = new TH1D(Form("h_llh_i_%s",GetName()),"log(p/pmax) interval",10000,-50,50);
  
  double lhmax   = fProb[fIPMax];
  double llh_bin = fHist.fLlh->GetXaxis()->GetBinWidth(0);
  
  for (int i=0; i<MaxNx; i++) {
    fHist.fProb->SetBinContent(i+1,fProb[i]);
    fHist.fProb->SetBinError  (i+1,0);

    double llh = -log(fProb[i])+log(lhmax);
    if (i < (fMuB+fMuS)) llh = -llh;
                                        // a hack: handle case of integer mean to avoid two entries in the same bin
                                        // try to move by one bin...
    if (i <  fIPMax) llh = llh-llh_bin;
      
    fHist.fLlh->Fill(llh,fProb[i]);

    if ((i >= fIxMin) and (i <= fIxMax)) {
      fHist.fInterval->SetBinContent(i+1,fProb[i]);
      fHist.fInterval->SetBinError  (i+1,0);

      fHist.fLlhInterval->Fill(llh,fProb[i]);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
void TBeltLH::make_belt_hist() {

  if (fHist.fBeltLo) {
    delete fHist.fBeltLo;
    delete fHist.fBeltHi;
    delete fHist.fBeltNO;
    delete fHist.fBeltUL;
  }
  
  fHist.fBeltLo   = new TH1D(Form("h_belt_lo_%s",GetName()),
                             Form("TBeltLH LO MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltHi   = new TH1D(Form("h_belt_hi_%s",GetName()),
                             Form("TBeltLH HI MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltNO   = new TH1D(Form("h_belt_no_%s",GetName()),
                             Form("TBeltLH NO MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltNO1   = new TH1D(Form("h_belt_no1_%s",GetName()),
                             Form("TBeltLH NO1 MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltUL   = new TH1D(Form("h_belt_ul_%s",GetName()),
                             Form("TBeltLH UL MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  for (int ix=0; ix<fBelt.fLlhNPoints; ix++) {
    fHist.fBeltLo ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix  ]);
    fHist.fBeltHi ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+1]);
    fHist.fBeltNO ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+2]);
    fHist.fBeltUL ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+3]);
    fHist.fBeltNO1->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+4]);
  }

  fHist.fBeltLo->GetXaxis()->SetTitle("#mu_{S}");
  fHist.fBeltLo->GetYaxis()->SetTitle("LLH");
  
  fHist.fBeltLo->SetFillStyle(fBelt.fFillStyle);
  fHist.fBeltLo->SetFillColor(fBelt.fFillColor);
  fHist.fBeltLo->SetLineColor(fBelt.fFillColor);

  fHist.fBeltHi->SetFillStyle(fBelt.fFillStyle);
  fHist.fBeltHi->SetFillColor(fBelt.fFillColor);
  fHist.fBeltHi->SetLineColor(fBelt.fFillColor);

  fHist.fBeltUL->SetFillStyle(3002);
  fHist.fBeltUL->SetFillColor(fBelt.fFillColor);
  fHist.fBeltUL->SetLineColor(fBelt.fFillColor);

  fHist.fBeltNO->SetMarkerStyle(20);
}


//-----------------------------------------------------------------------------
// assume S in 
//-----------------------------------------------------------------------------
int TBeltLH::test_coverage(double MuB, double SMin, double SMax, int NPoints) {
  
  int rc(0);
  rc = construct_belt(MuB,0,35,35001,-1,nullptr);
  if (rc < 0) return rc;

  float x[NPoints+2], y[NPoints+2];

  x[0] = SMin;
  y[0] = 0;
  double dy = (SMax-SMin)/(NPoints-1);
  for (int i=0; i<NPoints; i++) {
    double s = SMin + i*dy;
    if (s == 0) s = 1.e-10;          // deal with a numerical issue around zero
                                        // now generate pseudoexperiments
    int nmissed = 0;
    for (int k=0; k<fNExp; k++) {
      int nobs = fRn.Poisson(s+MuB);
                                        // assume SMin and SMax define a confidence interval
      double mus = nobs;
                                        // determine SMin and SMax;
      double smin = fBelt.fSign[nobs][0];
      double smax = fBelt.fSign[nobs][1];

      if ((s < smin) or (s > smax)) {
        nmissed += 1;
                                        // print only missed ones
        if (fDebug.fTestCoverage == 2) {
          printf("k, s, MuB, nobs, mus, smin, smax, nmissed : %10i %10.5f %10.3f %3i %10.3f %10.3f %10.3f %10i\n",
                 k, s, MuB, nobs, mus, smin, smax, nmissed);
        }
      }
      if (fDebug.fTestCoverage == 1) {
        printf("k, s, MuB, nobs, mus, smin, smax, nmissed : %10i %10.5f %10.3f %3i %10.3f %10.3f %10.3f %10i\n",
               k, s, MuB, nobs, mus, smin, smax, nmissed);
      }
    }
    float prob = 1.-nmissed/float(fNExp);
    x[i+1] = s;
    y[i+1] = prob;
  }

  x[NPoints+1] = SMax;
  y[NPoints+1] = 0;
  
  if (fHist.fCoverage) delete fHist.fCoverage;
  fHist.fCoverage = new TGraph(NPoints+2,x,y);
  fHist.fCoverage->SetTitle(Form("Crow-Gardner coverage test #mu_{B} = %10.3f",MuB));
  
  return rc;
}
}

