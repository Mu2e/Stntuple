///////////////////////////////////////////////////////////////////////////////
// my implementation of the Feldman-Cousins algorithm
// fDebug.fTestCoverage   = 1: print all
// fDebug.fTestCoverage   = 2: print missed
// fDebug.fConstructBelt >= 2: print within the range
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/stat/TFeldmanCousins.hh"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TROOT.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/DistFunc.h"

ClassImp(stntuple::TFeldmanCousins)
//-----------------------------------------------------------------------------
// CL > 0: CL has the meaning of a probability - 0.9, 0.95 .. .0.99 etc
//    < 0: CL is the "discovery probability", probability corresponding
//         to 5 gaussian sigma level, 
//-----------------------------------------------------------------------------
namespace stntuple {
TFeldmanCousins::TFeldmanCousins(const char* Name, double CL, int DebugLevel):
  TNamed(Name,Name),
  fRn()
{
  SetCL(CL);

  fType                     = 0;        // unbiased
  fDebug.fAll               = DebugLevel;
  fDebug.fConstructBelt     = 0;
  fDebug.fConstructInterval = 0;
  fDebug.fUpperLimit        = 0;
  fDebug.fTestCoverage      = 0;
  fDebug.fMuMin             = 0.;       // nothind printed by default
  fDebug.fMuMax             = -1.;
//-----------------------------------------------------------------------------
// calculate factorials - just once
//-----------------------------------------------------------------------------
  fFactorial[0] = 1;
  for (int i=1; i<MaxNx; i++) {
    fFactorial[i] = fFactorial[i-1]*i;
  }
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  fHist.fBgProb    = new TH1D(Form("h_bg_prob_%s"   ,GetName()),"h_bg_prob"   ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fBsProb    = new TH1D(Form("h_bs_prob_%s"   ,GetName()),"h_bs_prob"   ,MaxNx,-0.5,MaxNx-0.5);
  // fHist.fCumBgProb = new TH1D(Form("h_cumbg_prob_%s",GetName()),"h_cumbg_prob",MaxNx,-0.5,MaxNx-0.5);
  // fHist.fCumBsProb = new TH1D(Form("h_cumbs_prob_%s",GetName()),"h_cumbs_prob",MaxNx,-0.5,MaxNx-0.5);

  fHist.fProb      = new TH1D(Form("h_prob_2D_%s"    ,GetName()),"h prob 2D"  ,MaxNx,-0.5,MaxNx-0.5);
  fHist.fBelt      = nullptr;
  fHist.fBeltLo    = nullptr;
  fHist.fBeltHi    = nullptr;
  fHist.fCoverage  = nullptr;

  fBelt.fFillColor = kBlue+2;
  fBelt.fFillStyle = 3004;
  
  fNExp            = 10000000;
					// make sure uninitialzied values don't make sense
  fBelt.fNy        = -1;
  fBelt.fSMin      = 1e6;
  fBelt.fSMax      = -1e6;
}

TFeldmanCousins::~TFeldmanCousins() {
  delete fHist.fBgProb;
  delete fHist.fBsProb;
  // delete fHist.fCumBgProb;
  // delete fHist.fCumBsProb;
  delete fHist.fProb;
  if (fHist.fBelt) delete fHist.fBelt;
}

void  TFeldmanCousins::SetCL(double CL) {
  double alpha = 1-TMath::Erf(5./sqrt(2));
  
  if (CL > 0) fCL = CL;
  else        fCL = 1-alpha/2; // always two-sided: 1-5.7330314e-07

  fLog1mCL        = log(1-fCL);
}

//-----------------------------------------------------------------------------
void TFeldmanCousins::InitPoissonDist(double MuB, double MuS, double* Prob, int NObs) {

  double mean = MuB+MuS;
  Prob[0]     = TMath::Exp(-mean);

  if (NObs < 0) {
    for (int i=1; i<MaxNx; i++) {
      Prob[i] = Prob[i-1]*mean/i;
    }
  }
  else {
//-----------------------------------------------------------------------------
// NObs > 0 biases the probability distribution
//-----------------------------------------------------------------------------
    double pb[MaxNx];

    double pbn = 0;
    for (int k=0; k<=NObs; k++) {
      double pbb = TMath::Exp(-MuB)*pow(MuB,k)/fFactorial[k];
      double pbs = TMath::Exp(-MuS)*pow(MuS,NObs-k)/fFactorial[NObs-k];

      pb[k]      = pbb*pbs;
      pbn       += pb[k];
    }
    
    for (int i=0; i<MaxNx; i++) {
      if (i <= NObs) pb[i] = pb[i]/pbn;
      else           pb[i] = 0;
    }
					// 'i' - bin in the constrained Poisson distribution
    for (int i=0; i<MaxNx; i++) {
      double pi = 0;
					// an experiment observed N events, use pb[k]
      for (int k=0; k<=i; k++) {
	double ps = TMath::Exp(-MuS)*pow(MuS,i-k)/fFactorial[i-k];
	pi        = pi + pb[k]*ps;
      }
      Prob[i] = pi;
    }
  }
}

// //-----------------------------------------------------------------------------
// void TFeldmanCousins::Init(double Bgr, double Sig) {
//   //  fSigMean = SigMean;

//   fMuB = Bgr;
//   fMuS = Sig;
//                                         // init poisson doesn't redefine MuB and MuS
//   InitPoissonDist(fMuB     , fBgProb, -1);
//   InitPoissonDist(fMuB+fMuS, fBsProb, -1);
// }

//-----------------------------------------------------------------------------
  int TFeldmanCousins::ConstructInterval(double MuB, double MuS, int NObs) {
//-----------------------------------------------------------------------------
// output: [fIxMin,fIxMax] : a CL interval constructed using FC ordering for given 
//         MuB and MuS
//-----------------------------------------------------------------------------
  int rc(0);				// return code

  double prob[MaxNx], best_prob[MaxNx];

  if (fDebug.fAll > 10) {
    printf("TFeldmanCousins::ConstructInterval: MuB = %10.3f MuS = %10.3f\n",MuB,MuS);
  }
  
  fMuB = MuB;
  fMuS = MuS;
  fNObs = NObs;
                                        // init poisson doesn't redefine MuB and MuS
  InitPoissonDist(fMuB,   0, fBgProb, NObs);
  InitPoissonDist(fMuB,fMuS, fBsProb, NObs);
  
  for (int ix=0; ix<MaxNx; ix++) {
                                        // 'ix' : N(observed events)
    double sbest = ix-MuB;

    if (sbest <= 0) sbest = 0;
                                        // sb = ix, ALWAYS ! 
    double sb = sbest+MuB;

    fBestProb[ix] = TMath::Power(sb,ix)*TMath::Exp(-sb)/fFactorial[ix];
    fBestSig [ix] = sbest;

    if (NObs == -1) {
      fLhRatio [ix] = fBsProb[ix]/fBestProb[ix];
    }
    else if (NObs == -2) {
//-----------------------------------------------------------------------------
// measured NObs events, if Nobs < MuB, reduce P(best)
//-----------------------------------------------------------------------------
      fLhRatio [ix] = fBsProb[ix]/fBestProb[ix];
      double sbest  = ix-MuB;
      if (sbest < 0) fLhRatio [ix] *= exp(-sbest);
    }
    else {
                                        // biased
      
      InitPoissonDist(fMuB,fMuS , prob     , NObs);
      InitPoissonDist(fMuB,sbest, best_prob, NObs);
 
      fLhRatio [ix] = prob[ix]/best_prob[ix];
   }
  }
  
  if (fDebug.fAll > 10) {
    PrintData("BestProb" ,'d',fBestProb     ,18);
    PrintData("BestSig"  ,'d',fBestSig      ,18);
    PrintData("LhRatio"  ,'d',fLhRatio      ,18);
  }
//-----------------------------------------------------------------------------
// sort ranks
//-----------------------------------------------------------------------------
  double rmax;
    
  for (int i=0; i<MaxNx; i++) fRank[i] = i;

  for (int i1=0; i1<MaxNx-1; i1++) {
    rmax = fLhRatio[fRank[i1]];
    for (int i2=i1+1; i2<MaxNx; i2++) {
      double r2 = fLhRatio[fRank[i2]];
      if (r2 > rmax) {
	rmax      = r2;
	int i     = fRank[i1];
	fRank[i1] = fRank[i2];
	fRank[i2] = i;
      }
      if (fDebug.fConstructInterval > 0) {
        printf(" ----------- i1, i2: %3i %3i\n",i1,i2);
        PrintData("LhRatio" ,'d',fLhRatio  ,18);
        PrintData("Rank"    ,'i',fRank     ,18);
      }
    }
  }
  if (fDebug.fAll > 10) {
    PrintData("Rank   " ,'i',fRank     ,18);
  }
//-----------------------------------------------------------------------------
// build confidence interval corresponding to the probability fCL - 0.9 , 0.95, etc
//-----------------------------------------------------------------------------
  fIxMin       = MaxNx;
  fIxMax       = -1;

  int covered = 0;
  fProb       = 0;
  
  if (fDebug.fAll > 10) {
    printf(" ix fRank[ix] fBsProb[ind] fBestProb[ind] LhRatio[ind]   fProb     1-fProb  fIxMin fIxMax\n");
    printf(" --------------------------------------------------------------------------------------\n");
  }
  for (int ix=0; ix<MaxNx; ix++) {
    int ind = fRank[ix];
    fProb  += fBsProb[ind];
    if (ind < fIxMin) fIxMin = ind;
    if (ind > fIxMax) fIxMax = ind;

    if (fDebug.fAll > 10) {
      printf("%3i %8i   %10.3e     %10.3e  %10.3e  %10.3e %10.3e %5i %5i\n",
	     ix,ind,fBsProb[ind],fBestProb[ind],fLhRatio[ind],fProb,1-fProb,fIxMin,fIxMax);
    }
    
    if (fProb > fCL) {
      covered = 1;
      break;
    }
  }

  if (covered == 0) {
    if (fDebug.fAll > 0) {
      printf("TFeldmanCousins::ConstructInterval:TROUBLE: MuB, MuS: %12.5e %12.5e ",MuB,MuS);
      printf("prob = %12.5e fIxMin:%3i, fIxMax:%3i, 1-CL = %12.5e\n",fProb,fIxMin,fIxMax,1-fCL);
    }
    rc = -1;
  }
  else {
    if (fDebug.fAll > 0) {
      printf("TFeldmanCousins::ConstructInterval:success: MuB, MuS: %12.5e %12.5e ",MuB,MuS);
      printf("prob = %12.5e fIxMin:%3i, fIxMax:%3i, 1-CL = %12.5e\n",fProb,fIxMin,fIxMax,1-fCL);
    }
  }

  return rc;
}

//-----------------------------------------------------------------------------
int TFeldmanCousins::ConstructInterval(model_t* Model) {
//-----------------------------------------------------------------------------
// Model already comes with initialized probability distributions and cumulative probability
// distribution
//
// output: [fIxMin,fIxMax] : a CL interval constructed using FC ordering for given 
//         MuB and MuS
//-----------------------------------------------------------------------------
  int rc(0);				// return code

  double mub = Model->MuB();
  
  for (int ix=0; ix<MaxNx; ix++) {
    double sbest = ix-mub;

    if (sbest <= 0) sbest = 0;

    double sb = sbest+mub;

    fBestProb[ix] = TMath::Power(sb,ix)*TMath::Exp(-sb)/fFactorial[ix];
    fBestSig [ix] = sbest;
    fLhRatio [ix] = fBsProb[ix]/fBestProb[ix];
  }
//-----------------------------------------------------------------------------
// sort ranks
//-----------------------------------------------------------------------------
  double rmax;
    
  for (int i=0; i<MaxNx; i++) fRank[i] = i;

  for (int i1=0; i1<MaxNx-1; i1++) {
    rmax = fLhRatio[fRank[i1]];
    for (int i2=i1+1; i2<MaxNx; i2++) {
      double r2 = fLhRatio[fRank[i2]];
      if (r2 > rmax) {
	rmax      = r2;
	int i     = fRank[i1];
	fRank[i1] = fRank[i2];
	fRank[i2] = i;
      }
      if (fDebug.fConstructInterval > 0) {
        printf(" ----------- i1, i2: %3i %3i\n",i1,i2);
        PrintData("LhRatio" ,'d',fLhRatio  ,18);
        PrintData("Rank"    ,'i',fRank     ,18);
      }
    }

    if (fDebug.fAll > 0) {
      PrintData("Rank   " ,'i',fRank     ,18);
    }
  }
//-----------------------------------------------------------------------------
// build confidence interval corresponding to the probability fCL - 0.9 , 0.95, etc
//-----------------------------------------------------------------------------
  fIxMin       = MaxNx;
  fIxMax       = -1;

  int covered = 0;
  fProb       = 0;
  
  for (int ix=0; ix<MaxNx; ix++) {
    int ind = fRank[ix];
    fProb  += fBestProb[ind];
    if (ind < fIxMin) fIxMin = ind;
    if (ind > fIxMax) fIxMax = ind;

    if (fDebug.fAll > 0) {
      printf(" ix ind=fRank[ix] fBsProb[ind] fProb, 1-fProb fIxMin fIxMax :%3i %3i %10.3e %10.3e %10.3e %3i %3i\n",
	     ix,ind,fBsProb[ind],fProb,1-fProb,fIxMin,fIxMax);
    }
    
    if (fProb > fCL) {
      covered = 1;
      break;
    }
  }

  if (covered == 0) {
    printf("trouble ! interval not covered : prob = %12.5e , 1-CL = %12.5e\n",fProb,1-fCL);
  }

  return rc;
}

//-----------------------------------------------------------------------------
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
//-----------------------------------------------------------------------------
  int TFeldmanCousins::ConstructBelt(double MuB, double SMin, double SMax, int NPoints, int NObs) {

    if ((fBelt.fBgr == MuB) and (fBelt.fNy == NPoints) and (fBelt.fSMin == SMin) and (fBelt.fSMax == SMax)
        and (fBelt.fNObs == NObs)) return 0;

  fBelt.fBgr  = MuB;
  fBelt.fSMin = SMin;
  fBelt.fSMax = SMax;
  fBelt.fNy   = NPoints;
  fBelt.fNObs = NObs;
  
  fBelt.fDy   = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 1;

  for (int ix=0; ix<MaxNx; ix++) {
    fBelt.fSign[ix][0] = 1.e6;
    fBelt.fSign[ix][1] = -1;
  }

  for (int iy=0; iy<NPoints; iy++) {
    double mus = SMin+iy*fBelt.fDy;

    int rc     = ConstructInterval(MuB,mus,NObs);
    if (rc == 0) {
      for (int ix=fIxMin; ix<=fIxMax; ix++) {
        if (mus > fBelt.fSign[ix][1]) fBelt.fSign[ix][1] = mus+fBelt.fDy/2;
        if (mus < fBelt.fSign[ix][0]) fBelt.fSign[ix][0] = mus-fBelt.fDy/2;
        if (fDebug.fConstructBelt >= 2) {
          if ((mus >= fDebug.fMuMin) and (mus <= fDebug.fMuMax)) {
            printf("TFeldmanCousins::ConstructBelt: iy = %5i ix:%3i mus=%8.4f MuB=%5.3f IxMin:%3i IxMax:%3i fSign[ix][0]:%8.4f fSign[ix][1]:%8.4f\n",
                   iy,ix,mus,MuB,fIxMin,fIxMax,fBelt.fSign[ix][0],fBelt.fSign[ix][0]);
          }
        }
      }
    }
    else {
      if (fDebug.fConstructBelt > 0) {
	printf("TFeldmanCousins::ConstructBelt: ERROR: MuB=%10.4f mus=%10.4f IY = %3i interval not defined\n",
	       MuB, mus,iy);
      }
    }
  }

  for (int ix=0; ix<MaxNx; ix++) {
    if (fBelt.fSign[ix][0] == 1.e6) fBelt.fSign[ix][0] = 0;
    if (fBelt.fSign[ix][0] <  0   ) fBelt.fSign[ix][0] = 0;
    if (fBelt.fSign[ix][1] <  0   ) fBelt.fSign[ix][1] = 0;
  }

  if (fDebug.fConstructBelt > 0) {
    printf("       N(obs)      S(min)      S(max) \n");
    printf("------------------------------------- \n");
    for (int ix=0; ix<MaxNx; ix++) {
      printf("%10i %12.5f %12.5f\n",ix,fBelt.fSign[ix][0],fBelt.fSign[ix][1]);
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
//-----------------------------------------------------------------------------
void TFeldmanCousins::DiscoveryProb(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// construct FC CL confidence interval (covering fCL) assuming Signal=0
//-----------------------------------------------------------------------------
  ConstructInterval(MuB,0);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix]   = SMin+ix*step;
    double tot = MuB+MuS[ix];
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    long int ndisc = 0;			
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// definition of the discovery:
// rn > fIxMax, i.e  the  probability to observe'rn' is less than 1-fCL
//-----------------------------------------------------------------------------
      if (rn > fIxMax) ndisc ++;
    }
    Prob[ix] = double(ndisc)/double(fNExp);
  }
}

void TFeldmanCousins::DiscoveryProbMean(double MuB, double SMin, double SMax, int NPoints, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
// watch for 5 
// construct FC CL confidence interval (covering fCL) for MuS=0
//-----------------------------------------------------------------------------
  double mus;
  ConstructInterval(MuB,mus=0);

  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  for (int ix=0; ix<NPoints; ix++) {
    MuS[ix]    = SMin+ix*step;
    double tot = MuB+MuS[ix];
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    double sum  = 0;
    double sumn = 0;
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
      double p;
      if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
      else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);

      double      sig = ROOT::Math::gaussian_quantile_c(p,1);

      sum      += sig;
      sumn     += 1;
      if (fDebug.fAll > 0) {
	printf("i, tot, rn, p, sig, sum : %3i %10.4f %3i %15.8e %12.5e %12.5e\n",i,tot,rn, p,sig, sum);
      }
    }
    
    Prob[ix]  = sum/sumn; // ROOT::Math::gaussian_quantile(pp,1);

    if (fDebug.fAll > 0) {
      printf("ix, sum, sumn, MuS[ix], Prob[ix] : %3i %12.5e %10.3e %12.5e %12.5e\n",
	     ix,sum,sumn,MuS[ix],Prob[ix]);
    }
  }
}


void TFeldmanCousins::DiscoveryMedian(double MuB, double SMin, double SMax, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  *MuS  = -1;
  *Prob = -1;

  double NSIG(5.);
  double EPS(1.e-3);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  double ndisc = 0;
  
  double tot1 = MuB+SMin;
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(tot1);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
    if (sig > NSIG) ndisc++;
  }
    
  double p1 = ndisc/fNExp;

  printf("MuB = %10.5f Sig = %10.5f P = %10.5f\n",MuB,SMin,p1);


  ndisc = 0;
  double tot2 = MuB+SMax;
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(tot2);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
    if (sig > NSIG) ndisc++;
  }
    
  double p2 = ndisc/fNExp;
  printf("MuB = %10.5f Sig = %10.5f P = %10.5f\n",MuB,SMax,p2);

  if (((p1 < 0.5) and (p2 < 0.5)) or ((p1 > 0.5) and (p2 > 0.5))) {
    printf("eeeee\n");
    return;
  }

  double smin = SMin;
  double smax = SMax;

  double prob(-1);
  while (smax-smin > EPS) {
    double s   = (smin+smax)/2;
    double tot = s+MuB;

    ndisc = 0;
    for (int i=0; i<fNExp; i++) {
      int rn = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
      double p;
      if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
      else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);
      
      double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
      if (sig > NSIG) ndisc++;
    }

    prob = ndisc/fNExp;
    printf("MuB = %10.5f s   = %10.5f P = %10.5f\n",MuB,s,prob);

    if   (prob < 0.5) smin = s;
    else              smax = s;
  }

  *MuS  = (smin+smax)/2;
  *Prob = prob;

  printf("END: MuB = %10.5f s = %10.5f P = %10.5f\n",MuB,*MuS,*Prob);
}

//-----------------------------------------------------------------------------
void TFeldmanCousins::DiscoveryMedian(model_t* Model, double SMin, double SMax, double* MuS, double* Prob) {
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  *MuS  = -1;
  *Prob = -1;

  double sold = Model->SignalChannel()->Process()->Mean();

  double NSIG(5.  );
  double EPS(1.e-3);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  double ndisc = 0;

  double mub   = Model->GetBackgroundMean();
  
  double tot1 = mub+SMin;
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(tot1);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,mub);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,mub);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
    if (sig > NSIG) ndisc++;
  }
    
  double p1 = ndisc/fNExp;

  printf("MuB = %10.5f Sig = %10.5f P = %7.5f\n",mub,SMin,p1);


  ndisc = 0;
  double tot2 = mub+SMax;
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(tot2);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate to 'rn' or above
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,mub);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,mub);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
    if (sig > NSIG) ndisc++;
  }
    
  double p2 = ndisc/fNExp;
  printf("MuB = %10.5f Sig = %10.5f P = %7.5f\n",mub,SMax,p2); // 

  if (((p1 < 0.5) and (p2 < 0.5)) or ((p1 > 0.5) and (p2 > 0.5))) {
    printf("eeeee\n");
    return;
  }
//-----------------------------------------------------------------------------
// finally, do the job
//-----------------------------------------------------------------------------
  double smin = SMin;
  double smax = SMax;

  double mus, tot;
  int    nobs;

  double prob(-1);
  while (smax-smin > EPS) {
    double s   = (smin+smax)/2;
    Model->SignalChannel()->Process()->SetMean(s);

    ndisc = 0;
    for (int i=0; i<fNExp; i++) {
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
// next pseudoexperiment
//-----------------------------------------------------------------------------
      Model->InitParameters();
      
      mub  = Model->GetNullValue  ();
      mus  = Model->SignalChannel()->Process()->GetValue();
      tot  = mub+mus;
      
      nobs = fRn.Poisson(tot);

      double p;
      if (nobs > 0) p = ROOT::Math::poisson_cdf_c(nobs-1,mub);
      else          p = ROOT::Math::poisson_cdf_c(0     ,mub);
      
      double      sig = ROOT::Math::gaussian_quantile_c(p,1);
    
      if (sig > NSIG) ndisc++;
    }

    prob = ndisc/fNExp;
    printf("MuB = %10.5f s   = %10.5f P = %7.5f\n",mub,s,prob);

    if   (prob < 0.5) smin = s;
    else              smax = s;
  }
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
  Model->SignalChannel()->Process()->SetMean(sold);
  
  *MuS  = (smin+smax)/2;
  *Prob = prob;
  
  printf("END: MuB = %10.5f s = %10.5f P = %10.5f\n", Model->GetBackgroundMean(),*MuS,*Prob);
}


//-----------------------------------------------------------------------------
void TFeldmanCousins::DiscoveryProbMean(model_t* Model, double SMin, double SMax, int NPoints, double* MuS, double* NSig) {
//-----------------------------------------------------------------------------
// in general, need to scan a range of signals, call this function multiple times
// watch for 5 
// construct FC CL confidence interval (covering fCL) for MuS=0
//-----------------------------------------------------------------------------
  double step = (NPoints > 1) ? (SMax-SMin)/(NPoints-1) : 0;

  parameter_t* signal_process = Model->SignalChannel()->Process();

  for (int ix=0; ix<NPoints; ix++) {
					// scan a range of signal strengths
    MuS[ix]    = SMin+ix*step;
    signal_process->SetMean(MuS[ix]);	// should be fixed
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
    double sum  = 0;
    for (long int i=0; i<fNExp; i++) {
					// next pseudoexperiment: fluctuate nuisanse parameters
      Model->InitParameters();
					// fluctuated background mean for this pseudoexperiment
      
      double mub = Model->GetNullValue();
					// why do I need to construct an interval here? - probably, I don't
      //      ConstructInterval(Model);
					// signal strength, fluctuated in a correlated way with the background

      double mus = signal_process->GetValue();

      // model->GetValue() should do the same, but need to test

      // the total
      double tot = mub+mus;
      
      int    nobs = fRn.Poisson(tot);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate and result in an observatioin
// of 'nobs' or more events
//-----------------------------------------------------------------------------
      double p;
      if (nobs > 0) p = ROOT::Math::poisson_cdf_c(nobs-1,mub);
      else          p = ROOT::Math::poisson_cdf_c(0     ,mub);
//-----------------------------------------------------------------------------
// 'sig' is the deviation, in units of a 'gaussian sigma' corresponding
// to the probability 'p' of the background to fluctuate and give 'nobs' or more events
// such that p = Integral(sig,infinity)
//-----------------------------------------------------------------------------
      double sig = ROOT::Math::gaussian_quantile_c(p,1);

      sum      += sig;
      if (fDebug.fAll > 0) {
	printf("i, tot, nobs, p, sig, sum : %3li %10.4f %3i %15.8e %12.5e %12.5e\n",i,tot,nobs, p,sig, sum);
      }
    }
    
    NSig[ix]  = sum/fNExp; // ROOT::Math::gaussian_quantile(pp,1);

    if (fDebug.fAll > 0) {
      printf("ix, sum, fNExp, MuS[ix], NSig[ix] : %3i %12.5e %12li %12.5e %12.5e\n",
	     ix,sum,fNExp,MuS[ix],NSig[ix]);
    }
  }
}

//-----------------------------------------------------------------------------
void TFeldmanCousins::MakeProbHist() {

  if (fHist.fProb) {
    delete fHist.fProb;
    delete fHist.fBgProb;
    delete fHist.fBsProb;
  }
  
  fHist.fProb      = new TH1D(Form("h_prob_2D_%s"   ,GetName()),"h prob 2D",MaxNx,-0.5,MaxNx-0.5);
  fHist.fBgProb    = new TH1D(Form("h_bg_prob_%s"   ,GetName()),"h_bg_prob",MaxNx,-0.5,MaxNx-0.5);
  fHist.fBsProb    = new TH1D(Form("h_bs_prob_%s"   ,GetName()),"h_bs_prob",MaxNx,-0.5,MaxNx-0.5);
//-----------------------------------------------------------------------------
// [re]-initialize 1D histograms with the probabilities and integral probabilities
//-----------------------------------------------------------------------------
  for (int i=0; i<MaxNx; i++) {
    fHist.fProb  ->SetBinContent(i+1,fBsProb[i]);
    fHist.fBgProb->SetBinContent(i+1,fBgProb[i]);
    fHist.fBsProb->SetBinContent(i+1,fBsProb[i]);
  }
}
  
//-----------------------------------------------------------------------------
void TFeldmanCousins::MakeBeltHist() {
  if (fHist.fBelt) {
    delete fHist.fBelt;
    delete fHist.fBeltLo;
    delete fHist.fBeltHi;
  }
  
  fHist.fBeltLo   = new TH1D(Form("h_belt_lo_%s",GetName()),
                             Form("Feldman-Cousins belt lo MuB = %10.3f\n",fMuB),
                             MaxNx,-0.5,MaxNx-0.5);

  fHist.fBeltHi   = new TH1D(Form("h_belt_hi_%s",GetName()),
                             Form("Feldman-Cousins belt    MuB = %10.3f\n",fMuB),
                             MaxNx,-0.5,MaxNx-0.5);

  for (int ix=0; ix<MaxNx; ix++) {
    double dx = fBelt.fSign[ix][1]-fBelt.fSign[ix][0];
    if (fBelt.fSign[ix][1] > 0) {    
      fHist.fBeltLo->SetBinContent(ix+1,fBelt.fSign[ix][0]);
      fHist.fBeltHi->SetBinContent(ix+1,dx);
    }
  }

  fHist.fBeltLo->SetLineColor(fBelt.fFillColor);

  fHist.fBeltHi->SetFillStyle(fBelt.fFillStyle);
  fHist.fBeltHi->SetFillColor(fBelt.fFillColor);
  fHist.fBeltHi->SetLineColor(fBelt.fFillColor);

  fHist.fBelt = new THStack(Form("hs_%s",GetName()),fHist.fBeltHi->GetTitle());
  fHist.fBelt->Add(fHist.fBeltLo);
  fHist.fBelt->Add(fHist.fBeltHi);
}

//-----------------------------------------------------------------------------
void TFeldmanCousins::PlotDiscoveryProbMean(double MuB, double MuS) {
  // double mus;
  // ConstructInterval(MuB,mus=0);

  TH1F* h_sigm = new TH1F("h_sigm","h1 sigm"     ,2000,    0,20);
  TH1F* h_prob = new TH1F("h_prob","h1 log(prob)",2000, - 99, 1);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  if (fDebug.fAll > 0) {
    printf("  i  rn cum_prob[rn]        p          sig          sum       \n");
    printf("--------------------------------------------------------------\n");
  }
  for (int i=0; i<fNExp; i++) {
    int rn = fRn.Poisson(MuB+MuS);
    //    int rn = fRn.Poisson(MuS);
//-----------------------------------------------------------------------------
// define probability for the background to fluctuate above rn
//-----------------------------------------------------------------------------
    double p;
    if (rn > 0) p   = ROOT::Math::poisson_cdf_c(rn-1,MuB);
    else        p   = ROOT::Math::poisson_cdf_c(0   ,MuB);
    
    double      sig = ROOT::Math::gaussian_quantile_c(p,1);
      
    h_sigm->Fill(sig);
    h_prob->Fill(log(p));
  }

  const char* name = "c_fc_PlotDiscoveryProbMean";
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas("c_fc_PlotDiscoveryProbMean","c",1500,500);
    c->Divide(2,1);
  }
  c->cd(1);
  h_sigm->Draw();
  c->cd(2);
  h_prob->Draw();
}

//-----------------------------------------------------------------------------
void TFeldmanCousins::PrintData(const char* Title, char DataType, void* Data, int MaxInd) {

  int k      = 0;
  int findex = 0;
  int n_per_line(20);

  double*   pd = (double*) Data;
  int*      pi = (int*   ) Data;
  
  printf("-- %-10s: ",Title);
  for (int i=0; i<MaxInd; i++) {
    if (k == n_per_line) {
      printf("\n%03i",n_per_line*findex);
      k = 0;
    }
    if      (DataType == 'd') printf(" %9.2e",pd[i]);
    else if (DataType == 'i') printf(" %9i"  ,pi[i]);
    k++;
  }

  if (k > 0) printf("\n");
}

void TFeldmanCousins::PrintProbs(int N) {
//-----------------------------------------------------------------------------
// print results
//-----------------------------------------------------------------------------
  printf(" <bgr> = %10.4f <sig> = %10.4f\n"    ,fMuB,fMuS);
  printf(" IMin, IMax, Prob = %3i %3i %12.5f\n",fIxMin,fIxMax,fProb);

  PrintData("BgProb"  ,'d',fBgProb   ,N);
  PrintData("BestSig" ,'d',fBestSig  ,N);
  PrintData("BsProb"  ,'d',fBsProb   ,N);
  //  PrintData("CuBsProb",'d',fCumBsProb,N);
  PrintData("BestProb",'d',fBestProb ,N);
  PrintData("LhRatio" ,'d',fLhRatio  ,N);
  PrintData("Rank   " ,'i',fRank     ,N);
}

//-----------------------------------------------------------------------------
int TFeldmanCousins::SolveFor(double Val, const double* X, const double* Y, int NPoints, double* XVal) {
//-----------------------------------------------------------------------------
// fiven NPoints (X,Y) representing a smooth curve f(X), find XVal,
// such that f(XVal) = Val
//-----------------------------------------------------------------------------
  int i0 = -1;
  for (int i=0; i< NPoints; i++) {
    if (fDebug.fAll > 0) {
      printf("i, X[i], Y[i] : %2i %10.3e %12.5e\n",i,X[i],Y[i]);
    }
    if (Y[i] < Val) continue;
    i0 = i-1;
    break;
  }

  if (i0 == -1) {
    printf("TFeldmanCousins::solve_for in trouble: wrong range, bail out\n");
    return -1;
  }
  
  int i1 = i0+1;
  int i2 = i0+2;
  
  TMatrixD m(3,3);

  m(0,0) = X[i0]*X[i0]; m(0,1) = X[i0]; m(0,2) = 1;
  m(1,0) = X[i1]*X[i1]; m(1,1) = X[i1]; m(1,2) = 1;
  m(2,0) = X[i2]*X[i2]; m(2,1) = X[i2]; m(2,2) = 1;

  TMatrixD minv(TMatrixD::kInverted,m);

  TVectorD vy(3), vp(3);
  vy(0) = Y[i0];
  vy(1) = Y[i1];
  vy(2) = Y[i2];

  vp = minv*vy;

  double a = vp(0);
  double b = vp(1);
  double c = vp(2);
//-----------------------------------------------------------------------------
// finally, the solution for XVal
//-----------------------------------------------------------------------------
  if (a < 0) *XVal = -b/(2*a) - sqrt(b*b/(4*a*a)-(c-Val)/a);
  else       *XVal = -b/(2*a) + sqrt(b*b/(4*a*a)-(c-Val)/a);

  if (fDebug.fAll > 0) {
    printf ("a,b,c,x = %12.5e %12.5e %12.5e %12.5e \n",a,b,c,*XVal);
  }
  return 0;
}


//-----------------------------------------------------------------------------
  int TFeldmanCousins::TestCoverage(double MuB, double SMin, double SMax, int NPoints, int NObs) {
  
  int rc(0);
  rc = ConstructBelt(MuB,0,35,35001,NObs);

  if (rc < 0) return rc;

  float x[NPoints+2], y[NPoints+2];

  x[0] = SMin;
  y[0] = 0;

  double dy = (SMax-SMin)/(NPoints-1);
  
  for (int i=0; i<NPoints; i++) {
    double mus = SMin + i*dy;
    if (mus == 0) mus = 1.e-10;             // deal with a numerical issue around zero
                                        // now generate pseudoexperiments
    int nmissed = 0;
    for (int k=0; k<fNExp; k++) {
      int nobs = fRn.Poisson(mus+MuB);
                                        // determine SMin and SMax;
      double smin = fBelt.fSign[nobs][0];
      double smax = fBelt.fSign[nobs][1];

      if ((mus < smin) or (mus > smax)) {
        nmissed += 1;
                                        // print only missed ones
        if (fDebug.fTestCoverage == 2) {
          printf("k, mus, MuB, nobs, smin, smax, nmissed : %10i %10.5f %10.3f %3i %10.3f %10.3f %10i\n",
                 k, mus, MuB, nobs, smin, smax, nmissed);
        }
      }

      if (fDebug.fTestCoverage == 1) {
        printf("k, mus, MuB, nobs, smin, smax, nmissed : %10i %10.5f %10.3f %3i %10.3f %10.3f %10i\n",
               k, mus, MuB, nobs, smin, smax, nmissed);
      }
    }
    float prob = 1.- float(nmissed)/float(fNExp);
    x[i+1] = mus;
    y[i+1] = prob;
  }

  x[NPoints+1] = SMax;
  y[NPoints+1] = 0;

  if (fHist.fCoverage) delete fHist.fCoverage;
  fHist.fCoverage = new TGraph(NPoints+2,x,y);
  fHist.fCoverage->SetTitle(Form("coverage, MuB = %10.3f\n",MuB));
  
  return rc;
}

void TFeldmanCousins::UpperLimit(double MuB, double SMin, double SMax, int NPoints, double* S, double* Prob) {
//-----------------------------------------------------------------------------
// construct FC CL confidence belt
// constructing the FC belt assumes dividing the signal interval into (NPoints-1) steps,
// use 1000 points, assume signal in the range [0,10]
// assume MaxNy = 100+1, or similar, so MaxNy-1 is a multiple of 10
//-----------------------------------------------------------------------------
  ConstructBelt(MuB,0,10,1001);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  int nsteps = (NPoints > 1) ? NPoints -1 : 1;
  double step = (SMax-SMin)/nsteps;
  
//-----------------------------------------------------------------------------
// definition of exclusion: at least in (1-CL) of all cases, the value of S is above the belt
//-----------------------------------------------------------------------------
  for (int is=0; is<NPoints; is++) {
    S[is] = SMin+is*step;
    long int nexcl = 0;			
    for (int i=0; i<fNExp; i++) {
      int nobs = fRn.Poisson(MuB);
      if (S[is] > fBelt.fSign[nobs][1]) nexcl++;
      if (fDebug.fUpperLimit > 0) {
	printf("[TFeldmanCousins::UpperLimit] MuB, nobs S[is] fBelt.fSign[nobs][1] : %12.5e %3i %12.5e %12.5e\n",
	       MuB, nobs, S[is], fBelt.fSign[nobs][1]);
      }
    }
    Prob[is] = double(nexcl)/double(fNExp);
  }
}

int TFeldmanCousins::UpperLimit(double MuB, double SMin, double SMax, double* S, double* P) {
//-----------------------------------------------------------------------------
// construct FC CL confidence belt
// constructing the FC belt assumes dividing the signal interval into (NPoints-1) steps,
// use 1000 points, assume signal in the range [0,10]
// assume MaxNy = 100+1, or similar, so MaxNy-1 is a multiple of 10
//-----------------------------------------------------------------------------
  *S = -1;
  *P = -1;
  
  ConstructBelt(MuB,0,10,1001);
//-----------------------------------------------------------------------------
// ndisc: number of "discovery experiments", pseudoexperiments in which NULL
// hypothesis is excluded at (1-fCL) level
//-----------------------------------------------------------------------------
  long int nexcl = 0;
  for (int i=0; i<fNExp; i++) {
    int nobs = fRn.Poisson(MuB);
    if (SMin > fBelt.fSign[nobs][1]) nexcl++;
  }
  double p1 = double(nexcl)/double(fNExp);

  nexcl = 0;
  for (int i=0; i<fNExp; i++) {
    int nobs = fRn.Poisson(MuB);
    if (SMax > fBelt.fSign[nobs][1]) nexcl++;
  }
  double p2 = double(nexcl)/double(fNExp);

  if ((p1 > 0.5) or  (p2 < 0.5)) {
    // the point of interest is whole interval is beyond the interval
    printf("TFeldmanCousins::UpperLimit: CL beyond the interval: p(SMin), p(SMax): %10.5f %10.5f\n",p1,p2);
    return -1;
  }
					// can't do anything smarter as the function has breaks
  double p    (-1.);
  double EPS  (1.e-3);

  double smax = SMax;
  double smin = SMin;

  double s(-1);
  while ((p < 0.5) or (smax-smin > EPS)) {
    s = (smin+smax)/2;                 // calculate probability of S being excluded

    nexcl = 0;
    for (int i=0; i<fNExp; i++) {
      int nobs = fRn.Poisson(MuB);
      if (s > fBelt.fSign[nobs][1]) nexcl++;
    }

    p = double(nexcl)/double(fNExp);

    if (p > 0.5) smax = s;
    else         smin = s;

    printf(" s, smin, smax, p : %10.5f %10.5f %10.5f %10.5f\n",s,smin,smax,p);
  }
					// done
  *P = p;
  *S = s;

  return 0;
}


int TFeldmanCousins::UpperLimit(model_t* Model, double SMin, double SMax, double* S, double* P) {
//-----------------------------------------------------------------------------
// normally, used with CL = 90% or 95%
//-----------------------------------------------------------------------------
  *S = -1;
  *P = -1;
  
  double mub = Model->GetBackgroundMean();

  ConstructBelt(mub,0,20,401);
//-----------------------------------------------------------------------------
// nexcl: number of pseudoexperiments in which the NULL hypothesis is excluded at (1-fCL) level
// fCL fraction of pseudoexperiments excluded the hypothesis
//-----------------------------------------------------------------------------
  long int nexcl(0);
  for (int i=0; i<fNExp; i++) {
    int nobs = fRn.Poisson(mub);
    if (SMin > fBelt.fSign[nobs][1]) nexcl++;
  }
  double p1 = double(nexcl)/double(fNExp);

  nexcl = 0;
  for (int i=0; i<fNExp; i++) {
    int nobs = fRn.Poisson(mub);
    if (SMax > fBelt.fSign[nobs][1]) nexcl++;
  }
  double p2 = double(nexcl)/double(fNExp);

  printf("TFeldmanCousins::UpperLimit: p(SMin), p(SMax): %10.5f %10.5f\n",p1,p2);
  if ((p1 > 0.5) or  (p2 < 0.5)) {
    // the point of interest is whole interval is beyond the interval
    printf("TFeldmanCousins::UpperLimit: ERROR: beyond the interval, BAIL OUT\n");
    return -1;
  }
					// can't do anything smarter as the function has breaks
  double p    (-1.  );
  double EPS  (1.e-3);

  double smax = SMax;
  double smin = SMin;

  double s(-1);
  while ((p < 0.5) or (smax-smin > EPS)) {
    s = (smin+smax)/2;                 // calculate probability of S being excluded

    nexcl = 0;
    for (int i=0; i<fNExp; i++) {
      Model->InitParameters();		// next event

      mub = Model->GetNullValue();
					//

      if (fDebug.fUpperLimit > 0) {
	printf("i, mub : %10i %12.5f\n",i,mub);
      }
      
      ConstructBelt(mub,0,20,401);
      
      int nobs = fRn.Poisson(mub);
      if (s > fBelt.fSign[nobs][1]) nexcl++;
    }

    p = double(nexcl)/double(fNExp);

    if (p > 0.5) smax = s;
    else         smin = s;

    printf(" s, smin, smax, p : %10.5f %10.5f %10.5f %10.5f\n",s,smin,smax,p);
  }
					// done
  *P = p;
  *S = s;

  return 0;
}

  
};
