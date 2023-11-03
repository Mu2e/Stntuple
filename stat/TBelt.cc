///////////////////////////////////////////////////////////////////////////////
//
// fDebuglevel.fTestCoverage = 1: print all
// fDebuglevel.fTestCoverage = 2: print missed
///////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "Stntuple/stat/TBelt.hh"

ClassImp(stntuple::TBelt)

namespace stntuple {
//-----------------------------------------------------------------------------
TBelt::TBelt(const char* Name, double CL) : TNamed(Name,Name) {

  double alpha = 1-TMath::Erf(5./sqrt(2));
  
  if (CL > 0) fCL = CL;
  else        fCL = 1-alpha/2; // always one-sided: 5.7330314e-07/2 = 2.8665157e-07

  //  fType = Type;

  fIxMin = MaxNx;
  fIxMax = -1;
  
  fHist.fProb       = nullptr;
  fHist.fInterval   = nullptr;
  fHist.fBelt       = nullptr;
  fHist.fBeltLo     = nullptr;
  fHist.fBeltZr     = nullptr;
  fHist.fBeltHi     = nullptr;
  fHist.fCoverage   = nullptr;

  fDebug.fConstructBelt = 0;
  fDebug.fTestCoverage  = 0;
  fDebug.fUpperLimit    = 0;
                                        // no printout by default
  fDebug.fMuMin         = 1.;
  fDebug.fMuMax         = 0.;

  fNExp  = 100000;
  fNObs  = -1;

  fBelt.fFillColor      = kRed+2;
  fBelt.fFillStyle      = 3005;

  fBelt.fXMin           = 0;
  fBelt.fXMax           = -1;
  fBelt.fYMin           = 0;
  fBelt.fYMax           = -1;

  fBelt.fLlhInterval    = nullptr;
  fBelt.fLlhNPoints     = 0;
//-----------------------------------------------------------------------------
// calculate factorials, do that only once
// assume MaxNx to be large enough, so having N! values up to MaxNx-1 included is enough
//-----------------------------------------------------------------------------
  fFactorial[0] = 1; for (int i=1; i<MaxNx; i++) { fFactorial[i] = fFactorial[i-1]*i; }

  init_poisson_dist(0,1,fNObs);
}

//-----------------------------------------------------------------------------
// Zech version
//-----------------------------------------------------------------------------
int TBelt::init_truncated_poisson_dist(double MuB, int NObs, double* Prob) {

  double exp_mub  = TMath::Exp(-MuB);

  Prob[0]         = exp_mub;
  double sum_prob = Prob[0];
                                        // calculate normalization coefficient
  for (int k=1; k<=MaxNx; k++) {
    if (k <= NObs) {
      Prob[k]   = Prob[k-1]*MuB/k;
      sum_prob += Prob[k];
    }
    else {
      Prob[k] = 0;
    }
  }
                                        // do normalization
  for (int i=0; i<=NObs; i++) {
    Prob[i] = Prob[i]/sum_prob;
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
// Murat version
//-----------------------------------------------------------------------------
int TBelt::init_truncated_poisson_dist(double MuB, double MuS, int NObs, double* Prob) {

  printf("  TBelt::init_truncated_poisson_dist You don't want to be using this version!\n");
  return 0;
 
  double exp_mub  = TMath::Exp(-MuB);
  double exp_mus  = TMath::Exp(-MuS);

  double sum_prob = 0;
                                        // calculate normalization coefficients
  for (int k=0; k<=MaxNx; k++) {
    if (k <= NObs) {
      Prob[k]   = exp_mub*pow(MuB,k)/fFactorial[k] * (exp_mus*pow(MuS,NObs-k)/fFactorial[NObs-k]);
      sum_prob += Prob[k];
    }
    else {
      Prob[k] = 0;
    }
  }
                                        // do the normalization
  for (int i=0; i<=NObs; i++) {
    Prob[i] = Prob[i]/sum_prob;
  }
  
  return 0;
}


//-----------------------------------------------------------------------------
// mu = MuB+MuS
// 'NObs' is the number of observed events - it constrains the background fluctuations
// NObs<0 means no prior knowledge 
// array 'Prob' should have at least MaxNx elements
// declare NObs as double to be able to scan
// 'nature': given by MuB+Mus - sampled distribution
//-----------------------------------------------------------------------------
int TBelt::init_poisson_dist(double MuB, double MuS, int NObs) {
//-----------------------------------------------------------------------------
// the length of 'Prob' should be at least N
// CumProb[N]: probability, for a given Mean, to have rn <= N
//-----------------------------------------------------------------------------
  fMean  = MuB+MuS;
  fMuS   = MuS;
  fMuB   = MuB;
  fNObs  = NObs;
  
  if (NObs < 0) {
//-----------------------------------------------------------------------------
// standard Poisson distribution
//-----------------------------------------------------------------------------
    fProb[0] = TMath::Exp(-fMean);
    fIPMax   = 0;

    for (int i=1; i<MaxNx; i++) {
      fProb[i] = fProb[i-1]*fMean/i;
    }
    fIPMax = int(fMean);
  }
  else {
//-----------------------------------------------------------------------------
// background probability constrained by the measurement of N events (Zech'1989)
//-----------------------------------------------------------------------------
    double pb[MaxNx];

    // double pbn = 0; for (int k=0; k<=NObs; k++) { pbn += TMath::Exp(-MuB)*pow(MuB,k)/fFactorial[k]; }
    
    // for (int i=0; i<MaxNx; i++) {
    //   if (i <= NObs) pb[i] = TMath::Exp(-MuB)*pow(MuB,i)/fFactorial[i]/pbn;
    //   else           pb[i] = 0;
    // }

    init_truncated_poisson_dist(MuB,NObs,pb);
					// 'i' - bin in the constrained Poisson distribution
    double exp_mus = TMath::Exp(-MuS);

    double pmax = -1;
    fIPMax      = -1;
    for (int i=0; i<MaxNx; i++) {
					// experiment observed Nobs events, use pb[k]
      double pi = 0;
      for (int k=0; k<=i; k++) {
	double ps = exp_mus*pow(MuS,i-k)/fFactorial[i-k];
	pi        = pi + pb[k]*ps;
      }
      fProb[i] = pi;

      if (fProb[i] > pmax) {
	pmax   = fProb[i];
	fIPMax = i;
      }
    }
  }

  return 0;
}

  
//-----------------------------------------------------------------------------
// standard Poisson distribution
// the length of 'Prob' should be at least N
//-----------------------------------------------------------------------------
int TBelt::init_poisson_dist(double MuB, double* Prob) {

  fMean   = MuB;
  Prob[0] = TMath::Exp(-fMean);

  for (int i=1; i<MaxNx; i++) {
    Prob[i] = Prob[i-1]*fMean/i;
  }
  
  fIPMax = int(fMean);

  return 0;
}


  
//-----------------------------------------------------------------------------
// two parameters for uniformity of the interface
//-----------------------------------------------------------------------------
int TBelt::construct_interval(double MuB, double MuS) {
  printf("TBelt::construct_interval is undefined. BAIL OUT\n");
  return -1;
}

//-----------------------------------------------------------------------------
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
//-----------------------------------------------------------------------------
  int TBelt::construct_belt(double MuB, double SMin, double SMax, int NPoints, int NObs, double* P) {

  fMuB = MuB;
  
  fBelt.fSMin = SMin;
  fBelt.fSMax = SMax;
  fBelt.fDy   = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 1;

  for (int ix=0; ix<MaxNx; ix++) {
    fBelt.fSign[ix][0] = 1.e6;
    fBelt.fSign[ix][1] = -1;
  }

  for (int iy=0; iy<NPoints; iy++) {
    double mus = SMin+iy*fBelt.fDy;

    int rc     = construct_interval(MuB,mus);
    if (rc == 0) {
      for (int ix=fIxMin; ix<=fIxMax; ix++) {
        if (mus > fBelt.fSign[ix][1]) fBelt.fSign[ix][1] = mus+fBelt.fDy/2;
        if (mus < fBelt.fSign[ix][0]) fBelt.fSign[ix][0] = mus-fBelt.fDy/2;

        if (fDebug.fConstructBelt == 2) {
          if ((mus >= fDebug.fMuMin) and (mus <= fDebug.fMuMax)) {
            printf("TBelt::construct_belt: iy = %5i ix:%3i mus=%8.4f MuB=%5.3f IxMin:%3i IxMax:%3i fSign[ix][0]:%8.4f fSign[ix][1]:%8.4f\n",
                   iy,ix,mus,MuB,fIxMin,fIxMax,fBelt.fSign[ix][0],fBelt.fSign[ix][0]);
          }
        }
      }
    }
    else {
      if (fDebug.fConstructBelt > 0) {
	printf("TFeldmanCousinsB::ConstructBelt: ERROR: mus=%10.3f MuB=%10.3f iy = %3i interval not defined\n",
	       mus,MuB,iy);
      }
    }
  }

  for (int ix=0; ix<MaxNx; ix++) {
    if (fBelt.fSign[ix][0] == 1.e6) fBelt.fSign[ix][0] = 0;
    if (fBelt.fSign[ix][0] <  0   ) fBelt.fSign[ix][0] = 0;
    if (fBelt.fSign[ix][1] <  0   ) fBelt.fSign[ix][1] = 0;
  }

  if (fDebug.fConstructBelt > 0) {
    for (int ix=0; ix<MaxNx; ix++) {
      printf("%32s %3i %12.5f %12.5f\n","",ix,fBelt.fSign[ix][0],fBelt.fSign[ix][1]);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
int TBelt::make_prob_hist() {

  return 0;
}

//-----------------------------------------------------------------------------
void TBelt::make_belt_hist() {
  // if (fHist.fBelt) {
  //   delete fHist.fBelt;
  //   delete fHist.fBeltLo;
  //   delete fHist.fBeltHi;
  // }
  
  // fHist.fBeltLo   = new TH1D(Form("h_belt_lo_%s",GetName()),
  //                            Form("Crow-Gardner belt MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
  //                            MaxNx,-0.5,MaxNx-0.5);

  // fHist.fBeltZr   = new TH1D(Form("h_belt_mi_%s",GetName()),
  //                            Form("Crow-Gardner belt MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
  //                            MaxNx,-0.5,MaxNx-0.5);

  // fHist.fBeltHi   = new TH1D(Form("h_belt_hi_%s",GetName()),
  //                            Form("Crow-Gardner belt MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
  //                            MaxNx,-0.5,MaxNx-0.5);

  // for (int ix=0; ix<MaxNx; ix++) {
  //   double dx = fBelt.fSign[ix][1]-fBelt.fSign[ix][0];
  //   if (fBelt.fSign[ix][1] > 0) {    
  //     fHist.fBeltLo->SetBinContent(ix+1,fBelt.fSign[ix][0]);
  //     fHist.fBeltHi->SetBinContent(ix+1,dx);
  //   }
  // }

  // fHist.fBeltLo->GetXaxis()->SetTitle("N");
  // fHist.fBeltLo->GetYaxis()->SetTitle("#mu_{S}");
  
  // fHist.fBeltLo->SetLineColor(fBelt.fFillColor);
  // fHist.fBeltLo->SetLineColor(fBelt.fFillColor);

  // fHist.fBeltHi->SetFillStyle(fBelt.fFillStyle);
  // fHist.fBeltHi->SetFillColor(fBelt.fFillColor);
  // fHist.fBeltHi->SetLineColor(fBelt.fFillColor);

  // fHist.fBelt = new THStack(Form("hs_%s",GetName()),fHist.fBeltHi->GetTitle());
  // fHist.fBelt->Add(fHist.fBeltLo);
  // fHist.fBelt->Add(fHist.fBeltHi);
}


//-----------------------------------------------------------------------------
// assume S in 
//-----------------------------------------------------------------------------
int TBelt::test_coverage(double MuB, double SMin, double SMax, int NPoints) {
  
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
  fHist.fCoverage->SetTitle(Form("TBelt coverage test #mu_{B} = %10.3f",MuB));
  
  return rc;
}
}
