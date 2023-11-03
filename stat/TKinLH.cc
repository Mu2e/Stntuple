///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "TLine.h"

#include "Stntuple/stat/TKinLH.hh"

ClassImp(stntuple::TKinLH)

namespace stntuple {

int TKinLH::fDebug_QuickSort(0);

//-----------------------------------------------------------------------------
TKinLH::TKinLH(const char* Name, double CL, double PMin, double PMax, int Debug) : TBelt(Name, CL) {

  fInitialized              = 0;
  
  fDebug.fRun               = 0;
  fDebug.fConstructInterval = 0;
  fDebug.fConstructBelt     = 0;
  fDebug.fTestCoverage      = 0;
  fDebug.fMuMin             = 0;
  fDebug.fMuMax             = -1;
  pmin                      = PMin;
  pmax                      = PMax;
  fColor                    = kBlue+2;
  
  fHist.fProb       = nullptr;
  fHist.fLlh        = nullptr;
  fHist.fInterval   = nullptr;
  fHist.fBelt       = nullptr;
  fHist.fBeltLo     = nullptr;
  fHist.fBeltHi     = nullptr;
  fHist.fBeltProb   = nullptr;
  fHist.fBeltLhdt   = nullptr;
  fHist.fCoverage   = nullptr;

  fData.fNEvents    = 0;

  double tm3[] = {103., 104., 103., 104, 103., 104, 103., 104., 103., 104.};
  
  for (int i=0; i<10; i++) {
    fTestMom1[i] = 104.5;     // "signal-like"
    fTestMom2[i] = 102.0;     // "background-like"
    fTestMom3[i] = tm3[i];
  }

  init();
}

//-----------------------------------------------------------------------------
TKinLH::~TKinLH() {
  delete fSortData;
}

//-----------------------------------------------------------------------------
int TKinLH::generate() {
  return 0;
}

//-----------------------------------------------------------------------------
int TKinLH::init() {
//-----------------------------------------------------------------------------
// read in the CE momentum distribution and form the signal and the background
// probability distributions 'fHist.prob_sig' and 'fHist.prob_bgr'
// for the moment, consider uniform distribution for the background
//-----------------------------------------------------------------------------

  TH1F* h = gh1("/projects/mu2e/hist/su2020/su2020.cele0s61b1.su2020_track_ana.1010.hist",
                "su2020_TrackAna",
                "trk_2010/p");

  TString sig_name = Form("h_%s_prob_sig",GetName());
  TString bgr_name = Form("h_%s_prob_bgr",GetName());
    
  fHist.prob_sig = (TH1F*) h->Clone(sig_name.Data());
  fHist.prob_bgr = (TH1F*) h->Clone(bgr_name.Data());

  fHist.prob_bgr->Reset();

  double bin = fHist.prob_sig->GetBinWidth(1);

  int nx   = fHist.prob_sig->GetNbinsX();
  for (int i=0; i<nx; i++) {
    double x = fHist.prob_sig->GetBinCenter(i+1);
    if ((x < pmin) or (x > pmax)) {
      fHist.prob_sig->SetBinContent(i+1,0);
      fHist.prob_sig->SetBinError(i+1,0);
    }
    else {
      fHist.prob_bgr->SetBinContent(i+1,1);
      fHist.prob_bgr->SetBinError  (i+1,0);
    }
  }

  double sws = fHist.prob_sig->Integral();
  fHist.prob_sig->Scale(1./sws/bin);

  double swb = fHist.prob_bgr->Integral();
  fHist.prob_bgr->Scale(1./swb/bin);
//-----------------------------------------------------------------------------
// fSig is a fit of the signal probability distribution
// need that to avoid the binning effect
// in fact, fHist.prob_sig is not used for anything
//-----------------------------------------------------------------------------
  fSig = new TF1("f_sig",TKinLH::f_sig,pmin,pmax,1);
  fSig->SetParameter(0,1.);
  double scale = fSig->Integral(pmin,pmax);
  fSig->SetParameter(0,scale);

  fSig->SetNpx(30000);
  fSig->SetLineColor(kBlue+2);
//-----------------------------------------------------------------------------
// fBgr is a fit of the backgroundprobability distribution
// need that to avoid the binning effect
// in fact, fHist.prob_bgr is not used for anything
//-----------------------------------------------------------------------------
  fBgr = new TF1("f_bgr",TKinLH::f_bgr,pmin,pmax,1);
  fBgr->SetParameter(0,1.);
  scale = fBgr->Integral(pmin,pmax);
  fBgr->SetParameter(0,scale);

  fBgr->SetNpx(30000);
  fBgr->SetLineColor(kRed+2);
//-----------------------------------------------------------------------------
// book remaining histograms
// -------------------------
// 1. generated signal and background momenta (future event weights)
//-----------------------------------------------------------------------------
  TString name, title;
  
  title = Form("Generated P(bgr) name:%s pmin=%5.1f",GetName(),pmin);
  fHist.gen_pbgr = new TH1D("h_gen_pbgr",title,200,100,110);
    
  title = Form("Generated P(sig) name:%s pmin=%5.1f",GetName(),pmin);
  fHist.gen_psig = new TH1D("h_gen_psig",title,200,100,110);
  
  int     nbins_llh(20000), nbins_llhr(20000), nbins_llhrR(10000);
  
  double  llh_min  (-100), llh_max  (100);
  double  llhr_min (-100), llhr_max (100);

  for (int ntot=0; ntot<MaxNx; ntot++) {
                                        // nb varies from 0 to ntot
    
    fHist.fLogLhs [ntot]   = new TObjArray(ntot+1);
    fHist.fLogLhb [ntot]   = new TObjArray(ntot+1);
    fHist.fLogLhr [ntot]   = new TObjArray(ntot+1);
    fHist.fLogLhrR[ntot]   = new TObjArray(ntot+1);

    for (int nb=0; nb<=ntot; nb++) {
      TH1D* h;
      int ns = ntot-nb;
    
      title  = Form("Log(Lhs) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_llhs_%02i_%02i",nb,ns),title.Data(),nbins_llh,llh_min,llh_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhs[ntot]->Add(h);

      title  = Form("Log(Lhb) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_llhb_%02i_%02i",nb,ns),title.Data(),nbins_llh,llh_min,llh_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhb[ntot]->Add(h);

      title  = Form("Log(Lhr) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_llhr_%02i_%02i",nb,ns),title.Data(),nbins_llhr,llhr_min,llhr_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhr[ntot]->Add(h);

      title  = Form("Log(LhrR) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_llhrR_%02i_%02i",nb,ns),title.Data(),nbins_llhrR,0,1);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhrR[ntot]->Add(h);
    }
  }

  for (int ix=0; ix<MaxNx; ix++) {
    name  = Form("h_log_lhrR_1_n%02i",ix);
    title = Form("%s_log_lhrR_1_n%02i"  ,GetName(),ix);
    
    fHist.fLogLhrR_1[ix] = new TH1D(name,title,nbins_llhrR,0,50);

    name  = Form("h_log_lhrR_2_n%02i",ix);
    title = Form("%s_log_lhrR_2_n%02i"  ,GetName(),ix);
    
    fHist.fLogLhrR_2[ix] = new TH1D(name,title,nbins_llhrR,0,50);
  }

  name  = Form("h_sum_log_lhrR_1");
  title = Form("%s_sum_log_lhrR_1"  ,GetName());
  
  fHist.fSumLogLhrR_2 = new TH1D(name,title,nbins_llhrR,0,50);

  fSortData = new sdata[nbins_llhrR];
//-----------------------------------------------------------------------------
// all histograms booked, determine the LLHR range (for one event),
// scan the range with a small step,
// say, 10,000+1 steps, to include the interval ends
// this step is also model-dependent
//-----------------------------------------------------------------------------
  fMinLLHR    =  1e6;
  fMaxLLHR    = -1e6;
  int nsteps  = 10000;
  double step = (pmax-pmin)/nsteps;
    
  for (int i=0; i<nsteps+1; i++) {
    double p     = pmin+i*step;
                                        // define likelihoods only within the allowed momentum range
    double llhb  = log(lh_bgr(p));
    double llhs  = log(lh_sig(p));
    double llhr  = llhb-llhs;
    if (llhr < fMinLLHR) fMinLLHR = llhr;
    if (llhr > fMaxLLHR) fMaxLLHR = llhr;
  }

  return 0;
}

//-----------------------------------------------------------------------------
double TKinLH::bgr_mom() {
  double p = fBgr->GetRandom();
  return p;
}

//-----------------------------------------------------------------------------
double TKinLH::sig_mom() {
  double p = fSig->GetRandom();
  return p;
}

//-----------------------------------------------------------------------------
// normalized to the integral (not sum over the bins !)  = 1
//-----------------------------------------------------------------------------
double TKinLH::lh_bgr(double P) {
  double p = fBgr->Eval(P);
  return p;
}

//-----------------------------------------------------------------------------
// normalized to the integral (not sum over the bins !)  = 1
//-----------------------------------------------------------------------------
double TKinLH::lh_sig(double P) {
  double p = fSig->Eval(P);
  return p;
}


//-----------------------------------------------------------------------------
// normalized to the integral (not sum over the bins !)  = 1
// NObs : number of observed events, P - array of the track momenta
//-----------------------------------------------------------------------------
  double TKinLH::wt_data(double MuB, double MuS, int NObs, double* P, double* LlhrR) {
  
  double pb[MaxNx];

  int rc   = init_truncated_poisson_dist(MuB,NObs,pb);
  if (rc < 0) return rc;
//-----------------------------------------------------------------------------
// step 1: calculate reduced kinematic likelihood
//-----------------------------------------------------------------------------
  double tot_lhr = 1;
  if (NObs > 0) {
    for (int i=0; i<NObs; i++) {
      double p    = P[i];
      double lhr  = lh_bgr(p)/lh_sig(p);
      tot_lhr    *= lhr;
    }

    *LlhrR  = (log(tot_lhr)/NObs-fMinLLHR)/(fMaxLLHR-fMinLLHR);
  }
  else {
                                        // use the 1st bin of the uniform LogLhrR distribution for N=0
    *LlhrR = 0;
  }
//-----------------------------------------------------------------------------
// step 2: calculate probability of reduced kinematic likelihood
//         (assume the interval has just been constructed and fHist.fLogLhrR_1 filled)
//-----------------------------------------------------------------------------
  TH1D*  h   = fHist.fLogLhrR_1[NObs];
  int    bin = h->FindBin(NObs+*LlhrR);
                                        // 'wt' is the (binned) dP/d(llhrR) 
  double wt  = h->GetBinContent(bin);
//-----------------------------------------------------------------------------
// step 3: calculate global probability, for given MuB and MuS
//-----------------------------------------------------------------------------
  double exp_mus = TMath::Exp(-MuS);
  double p_nobs  = 0;
  for (int nb=0; nb<=NObs; nb++) {
    int ns    = NObs-nb;
    double ps = exp_mus*pow(MuS,ns)/fFactorial[ns];
    p_nobs   += pb[nb]*ps;
  }

  double prob = wt*p_nobs;

  return prob;
}


//-----------------------------------------------------------------------------
// P[0] - integral of the function over [pmin,pmax]
//-----------------------------------------------------------------------------
double TKinLH::f_bgr(double* X, double * P) {
  return 1./P[0];
};


double TKinLH::f_sig(double* X, double * P) {
  double f(0);

  int const nranges(4);

  double prange[nranges][2] = {
    102.00, 103.50,
    103.40, 104.20,
    104.15, 104.70,
    104.65, 105.10
  };

  double pmean(105);
  
  double par[nranges][4] = {
     2.60578e-02,  2.49680e-02,  8.57067e-03,  1.01647e-03,
     2.99588e-02,  2.66528e-02,  7.15730e-03,  4.32025e-04,
    -1.81856e-02, -1.31032e-01, -1.71238e-01, -7.00615e-02,
     5.95903e-04, -5.73379e-03,  3.24206e-02, -9.19974e-02
  } ;

  double p  = X[0];
  double dp = p-pmean;

  // 1. figure out the range, assume p>prange[0]

  if ((p < prange[0][0]) or (p > prange[nranges-1][1])) {
    printf(" momentum p = %10.3f is outside the range,  BAIL OUT\n",p);
    return 1;
  }

  int ir = -1;
  for (int i=0; i<nranges; i++) {
    if (p <= prange[i][1]) {
      ir = i;
      break;
    }
  }

  if ((ir < nranges-1) and (p >= prange[ir+1][0])) {
    // overlap region
    double f1 = par[ir  ][0] + par[ir  ][1]*dp + par[ir  ][2]*dp*dp + par[ir  ][3]*dp*dp*dp;
    double f2 = par[ir+1][0] + par[ir+1][1]*dp + par[ir+1][2]*dp*dp + par[ir+1][3]*dp*dp*dp;

    double pmin = prange[ir+1][0];
    double pmax = prange[ir  ][1];
    
    double dpp = pmax-pmin;
    
    f = (f2*(p-pmin) + f1*(pmax-p))/dpp;
  }
  else {
    // no overlap
    f = par[ir][0] + par[ir][1]*dp + par[ir][2]*dp*dp + par[ir][3]*dp*dp*dp;
  }

  return f/P[0];
}

//-----------------------------------------------------------------------------
int TKinLH::partition(sdata a[], int low, int high, double pivot) {

                                        // assume pivot = a[imax]
  int imax  = low;
  int i     = imax+1;

  while( i <= high){
    if (a[i].x <= pivot) {
      i++;
    }
    else{
      sdata temp;
      temp    = a[i];
      a[i]    = a[imax];
      a[imax] = temp;
      imax++;
      i++;
    }
  }
  
  return imax;
}

//-----------------------------------------------------------------------------
void  TKinLH::quickSort(sdata a[], int low, int high) {
  static int ncalls(0), l_low, l_high;
  
  if (low < high) {
    if (fDebug_QuickSort) {
      if (ncalls == 0) {
        l_low  = low;
        l_high = high;
      }
  
      printf("quickSort: ncalls: %i low, high: %5i %5i\n",ncalls, low,high);
      printf("TKinLH::quickSort 001: low:%2i high:%2i \n",low,high);
      int max_per_line(10), n_per_line(0);
    
      for (int i=l_low; i<=l_high; i++) {
        if (n_per_line == 0) printf("%5i:",i);
        printf(" (%2i,%9.2e)",a[i].bin,a[i].x);
        n_per_line++;
        if (n_per_line == max_per_line) {
          printf("\n");
          n_per_line = 0;
        }
      }
      if (n_per_line > 0) printf("\n");
      printf("TKinLH::quickSort 002: low:%2i high:%2i pivot:%9.2e\n",low,high,a[low].x);
    }
    
    int pos      = partition(a, low, high, a[low].x);
    
    ncalls++;
    
    if ((fDebug_QuickSort > 0) and (ncalls > fDebug_QuickSort)) return;
    
    quickSort(a,low  , pos-1);
    quickSort(a,pos+1, high );
  }
}

//-----------------------------------------------------------------------------
// three parameters, to maintain uniform interface
// Nobs is NOT USED, Zech was a mistake
// the step can be completed only after all log_lhrR histograms are filled
// assume the histograms are read in
//-----------------------------------------------------------------------------
int TKinLH::construct_interval(double MuB, double MuS) {

  double pb[MaxNx];

  int rc = init_poisson_dist(MuB,pb);
  
  if (fDebug.fConstructInterval) {
    printf("TKinLH::construct_interval 001:\n");
                                        // assume MaxNx % 10 = 0
    for (int i=0; i<MaxNx; i++) {
      printf(" %12.5e",pb[i]);
      if (((i+1) % 10) == 0) {
        printf("\n");
      }
    }
  }
//-----------------------------------------------------------------------------
// next: for given MuB and MuS, construct LogLhrR_N histograms
// LogLhrR_1: distribution in llhrR for a given ntot, summed over all nb with
//            proper weights (assuming known MuB abd MuS)
// LogLhrR_2: 
//-----------------------------------------------------------------------------
  double exp_mus = TMath::Exp(-MuS);
                                        // calculate P(NObs), assume MaxNx to be large enough

  TH1D*  h0      = (TH1D*) fHist.fLogLhrR[0]->At(0);
  int    nx      = h0->GetNbinsX();
//-----------------------------------------------------------------------------
// for given MuB and MuS, this loop coul dbe executed just once
//-----------------------------------------------------------------------------
  for (int nt=0; nt<MaxNx; nt++) {
    TObjArray* arR = fHist.fLogLhrR[nt];
    fHist.fLogLhrR_1[nt]->Reset();
    
    double p_nt = 0;
    for (int nb=0; nb<=nt; nb++) {
      int    ns = nt-nb;
      p_nt += pb[nb]*exp_mus*pow(MuS,ns)/fFactorial[ns];
    }
  
    if (fDebug.fConstructInterval) {
      printf("TKinLH::construct_interval 002: nt = %3i p_nobs = %12.5e\n",nt,p_nt);
    }
  
    for (int nb=0; nb<=nt; nb++) {
      int    ns = nt-nb;
      double ps = exp_mus*pow(MuS,ns)/fFactorial[ns];
      
                                        // this is the absolute normalization of the corresponding histogram
      double pns = pb[nb]*ps;
      
      TH1D* h     = (TH1D*) arR->At(nb);
                                        // kludge: make sure the hist is normalized, should've been already done
      double total = h->Integral();
      h->Scale(1./total);

      for (int ix=0; ix<nx; ix++) {
        double llhrR = h->GetBinCenter (ix+1);
//-----------------------------------------------------------------------------
// global llhR : can only assign observables - correct by 'p_nobs', not 'pns'
// the probability, however is defined by 'pns'
//-----------------------------------------------------------------------------
        double llhrG = llhrR+nt;
        double wt    = h->GetBinContent(ix+1)*pns;
                                        // logLhrR_N is normalized to the Poisson probability P(MuB,MuS,NObs)
                                        // just summing over all hists with the same NObs
        fHist.fLogLhrR_1[nt]->Fill(llhrG,wt);
      }
    }
  }
//-----------------------------------------------------------------------------
// now, create uniformly normalized distributions, 2-sided
//-----------------------------------------------------------------------------
  for (int nt=0; nt<MaxNx; nt++) {
    fHist.fLogLhrR_2[nt]->Reset();
//-----------------------------------------------------------------------------    
// h1 is supposed to be normalized to an integral (sum of contens of all bins) of 1
//-----------------------------------------------------------------------------
    TH1D*  h1 = fHist.fLogLhrR_1[nt];

    for (int ib=0; ib<nx; ib++) { 
      double wt    =  h1->GetBinContent(ib+1); // /pmax;
      double llhrG =  h1->GetBinCenter (ib+1);
      // if (nt >= (MuB+MuS)) llhrG = -llhrG;
      fHist.fLogLhrR_2[nt]->Fill(llhrG,wt);
    }
  }
//-----------------------------------------------------------------------------
// as a cross-check, fHist.fSumLogLhrR_2 should be normalized to unity
//-----------------------------------------------------------------------------
  TH1D* h_sum = fHist.fSumLogLhrR_2;
  nx          = h_sum->GetNbinsX();
  h_sum->Reset();
  
  for (int nt=0; nt<MaxNx; nt++) {
    h_sum->Add(fHist.fLogLhrR_2[nt]);
  }
                                        // set errors to 0
  for (int ib=0; ib<nx; ib++) {
    h_sum->SetBinError(ib+1,0);
  }
//-----------------------------------------------------------------------------
// last step: define the interval in the likelihod_ratio space. remember - it is two-sided
// everything starts from bin nx/2+1
// assume nb is an even number
//-----------------------------------------------------------------------------
  for (int ib=0; ib<nx; ib++) {
    fSortData[ib].bin = ib+1;
    fSortData[ib].x   = h_sum->GetBinContent(ib+1);
  }

  if (fDebug.fConstructInterval) {
    printf("TKinLH::construct_interval 003: before sorting, nx = %i\n",nx);
  }
  
  quickSort(fSortData,0,nx-1);

  if (fDebug.fConstructInterval) {
    printf("TKinLH::construct_interval 004: done sorting\n");
  }
                                        // defined the interval
  double sump         = 0;
  int    overcoverage = 0;
  double pmin         = 1.e6;
  
  for (int i=0; i<nx; i++) {
    int    bin = fSortData[i].bin;
    double p   = fSortData[i].x;
    sump       = sump+p;
    if (sump >= fCL) {
      if (overcoverage == 0) {
	if (i > 0) pmin = fSortData[i-1].x;
	if (p < pmin) {
					// the probability density is lower then in the previous bin
					// done, this bin needs to be included
	  fInterval.fLlhrMin = 0;
	  fInterval.fLlhrMax = h_sum->GetBinCenter(bin);  // interval bound - always positive
	  fInterval.fProbTot = sump;
	  fInterval.fPMax    = p;
	  fInterval.fIMax    = i+1;       // bins included into the region
	  break;
	} 
	else {
					// in the previous bin, the probability density was the same
					// this bin needs to be included , as the others with the same 
					// probability density, but we're overcovering
	  overcoverage = 1;
	}
      }
      else {
					// in the "overcoverage" mode continue till the first bin with 
					// lower prob
	if (p < fSortData[i-1].x) {
					// the probability density is lower then in the previous bin
					// done, this bin doesn't need to be included
	    fInterval.fLlhrMin = 0;
	    int    bin = fSortData[i-1].bin;
	    fInterval.fLlhrMax = h_sum->GetBinCenter(bin);  // interval bound - always positive
	    fInterval.fProbTot = sump-p;
	    fInterval.fPMax    = pmin;
	    fInterval.fIMax    = i;     // bin doesn't need to be included into the region
	    break;
	  } 
      }
                                        // otherwise, just redefined pmin
      pmin = p;
    }
  }
  
  if (fDebug.fConstructInterval) {
    printf("TKinLH::construct_interval 005: END\n");
  }
  
  return rc;
}
  
//-----------------------------------------------------------------------------
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
// NObs is only used to calculate the likelihood of the observation:
// 
//-----------------------------------------------------------------------------
  int TKinLH::construct_belt(double MuB, double SMin, double SMax, int NPoints, int NObs, double* P) {

  fMuB  = MuB;
  
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
    double mus   = SMin+i*fBelt.fDy;

    int rc       = construct_interval(MuB,mus);

    if (rc == 0) {
      fBelt.fLlhInterval[5*i  ] = fInterval.fLlhrMin;
      fBelt.fLlhInterval[5*i+1] = fInterval.fLlhrMax;
      fBelt.fLlhInterval[5*i+2] = fInterval.fProbTot;
      fBelt.fLlhInterval[5*i+3] = fInterval.fPMax;
//-----------------------------------------------------------------------------
// if data defined, also store the data line
//-----------------------------------------------------------------------------
      double wt = -1;
      if (P) {
        double llhrR;
        wt  = wt_data(MuB,mus,NObs,P,&llhrR);
      }
      fBelt.fLlhInterval[5*i+4] = -log(wt);
    }

    if (fDebug.fConstructBelt > 0) {
      for (int i=0; i<NPoints; i++) {
        printf("i,fBelt.fLlhInterval[5*i  ],fBelt.fLlhInterval[5*i+1]: %12.5f %12.5f %12.5f\n",
               fBelt.fLlhInterval[5*i  ],fBelt.fLlhInterval[5*i+1], fBelt.fLlhInterval[5*i+2]);
      }
    }
  }
  return 0;
}


//-----------------------------------------------------------------------------
void TKinLH::make_belt_hist() {

  if (fHist.fBelt) {
    delete fHist.fBelt;
    delete fHist.fBeltLo;
    delete fHist.fBeltHi;
    delete fHist.fBeltProb;
  }
  
  fHist.fBeltLo   = new TH1D("h_belt_lo",
                             Form("TBeltLH LO MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltHi   = new TH1D("h_belt_hi_%s",
                             Form("TBeltLH HI MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltProb = new TH1D("h_belt_prob",
                             Form("TBeltLH Prob MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltLhdt = new TH1D("h_belt_lhdt_%s",
                             Form("TBeltLH lhdt MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  for (int ix=0; ix<fBelt.fLlhNPoints; ix++) {
    double llhmax = -log(fBelt.fLlhInterval[5*ix+3]);
    double prob   = fBelt.fLlhInterval[5*ix+2];
    double lhdt   = fBelt.fLlhInterval[5*ix+4];

    fHist.fBeltLo  ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix  ]);
    fHist.fBeltHi  ->SetBinContent(ix+1,llhmax);
    fHist.fBeltProb->SetBinContent(ix+1,prob);
    fHist.fBeltLhdt->SetBinContent(ix+1,lhdt);
  }

  fHist.fBeltLo->GetXaxis()->SetTitle("#mu_{S}");
  fHist.fBeltLo->GetYaxis()->SetTitle("LLH");
  
  fHist.fBeltLo->SetLineColor(fBelt.fFillColor);

  fHist.fBeltHi->SetFillStyle(fBelt.fFillStyle);
  fHist.fBeltHi->SetFillColor(fBelt.fFillColor);
  fHist.fBeltHi->SetLineColor(fBelt.fFillColor);

  fHist.fBeltProb->SetMarkerStyle(20);
  fHist.fBeltLhdt->SetMarkerStyle(20);

  fHist.fBelt = new THStack(Form("hs_%s",GetName()),fHist.fBeltHi->GetTitle());
  fHist.fBelt->Add(fHist.fBeltLo);
  fHist.fBelt->Add(fHist.fBeltHi);
}

//-----------------------------------------------------------------------------
void TKinLH::plot_interval() {
  fHist.fSumLogLhrR_2->Draw("");

  TH1D* h1 = (TH1D*) fHist.fSumLogLhrR_2->Clone("h1");
  h1->Reset();

  int nx = fInterval.fIMax;
  
  for (int i=0; i<nx; i++) {
    int    bin  = fSortData[i].bin;
    double data = fSortData[i].x;

    // printf("bin, data: %5i %12.5e\n",bin,data);
    
    h1->SetBinContent(bin,data);
    h1->SetBinError  (bin,0);
  }

  h1->SetLineColor  (kRed+2);
  h1->SetMarkerColor(kRed+2);
  h1->SetFillStyle(3004);
  h1->SetFillColor(kRed+2);
  h1->Draw("sames");

  TLine* line = new TLine(-50,fInterval.fPMax,50,fInterval.fPMax);
  line->SetLineColor(kRed+2);
  line->Draw();
}


//-----------------------------------------------------------------------------
int TKinLH::read_hist(const char* Filename) {

  if (fInitialized) {
    printf("TKinLH::read_hist: re-initialization attempt, BAIL OUT\n");
    return -1;
  }
  
  TFile* f = TFile::Open(Filename);
  gROOT->cd();
                                        // current directory is gROOT
  
  for (int nobs=0; nobs<MaxNx; nobs++) {

    TObjArray* as  = fHist.fLogLhs [nobs];
    TObjArray* ab  = fHist.fLogLhb [nobs];
    TObjArray* ar  = fHist.fLogLhr [nobs];
    TObjArray* arR = fHist.fLogLhrR[nobs];

    for (int nb=0; nb<=nobs; nb++) {
      int ns = nobs-nb;
      (*as) [nb]->Delete();
      (*as) [nb] = (TH1D*) f->Get(Form("//%05i/h_llhs_%02i_%02i" ,nobs,nb,ns));
      (*ab) [nb]->Delete();
      (*ab) [nb] = (TH1D*) f->Get(Form("//%05i/h_llhb_%02i_%02i" ,nobs,nb,ns));
      (*ar) [nb]->Delete();
      (*ar) [nb] = (TH1D*) f->Get(Form("//%05i/h_llhr_%02i_%02i" ,nobs,nb,ns));
      (*arR)[nb]->Delete();
      (*arR)[nb] = (TH1D*) f->Get(Form("//%05i/h_llhrR_%02i_%02i",nobs,nb,ns));
    }

    // fHist.fLogLhrR_1[nobs]->Read(Form("h_llhrR1_%02i" ,nobs));
    // fHist.fLogLhrR_2[nobs]->Read(Form("h_llhrR2_%02i" ,nobs));
  }

  delete fHist.prob_sig;
  fHist.prob_sig = (TH1F*) f->Get("h_prob_sig");
  delete fHist.prob_bgr;
  fHist.prob_bgr = (TH1F*) f->Get("h_prob_bgr");

  // fHist.fSumLogLhrR_2->Read(Form("h_sum_llhrR2_%02i" ,nobs));

  // f->Close();
  
  // delete f;

  fInitialized = 1;
  
  return 0;
}
  

//-----------------------------------------------------------------------------
int TKinLH::run(int NObs, int NPe) {

  TObjArray* as  = fHist.fLogLhs [NObs];
  TObjArray* ab  = fHist.fLogLhb [NObs];
  TObjArray* ar  = fHist.fLogLhr [NObs];
  TObjArray* arR = fHist.fLogLhrR[NObs];

  for (int nb=0; nb<=NObs; nb++) {
    TH1D* h = (TH1D*) as->At(nb);
    h->Reset();
    h->SetLineColor  (fColor);
    h->SetMarkerColor(fColor);
    
    h = (TH1D*) ab->At(nb);
    h->Reset();
    h->SetLineColor  (fColor);
    h->SetMarkerColor(fColor);
      
    h = (TH1D*) ar->At(nb);
    h->Reset();
    h->SetLineColor  (fColor);
    h->SetMarkerColor(fColor);

    h = (TH1D*) arR->At(nb);
    h->Reset();
    h->SetLineColor  (fColor);
    h->SetMarkerColor(fColor);
  }

  fHist.gen_pbgr->Reset();
  fHist.gen_psig->Reset();

  if (NObs == 0) {
//-----------------------------------------------------------------------------
// special case - no kinematic info, have to make it up in a coherent way
// need a uniform distribution in a range defined by the probability distributions
// there is only one Lhr histogram, individual likelihoods are not defined 
//-----------------------------------------------------------------------------
    TH1D* h = (TH1D*) fHist.fLogLhr[0]->At(0);
    
    int nb = h->GetNbinsX();

    for (int ib=0; ib<nb; ib++) {
      double llhr = h->GetBinCenter(ib+1);
      if ((llhr >= fMinLLHR) and (llhr <= fMaxLLHR)) {
        h->SetBinContent(ib+1,1);
        h->SetBinError  (ib+1,0);
      }
    }
 
    TH1D* h1 = (TH1D*) fHist.fLogLhrR[0]->At(0);
    nb = h1->GetNbinsX();

    for (int ib=0; ib<nb; ib++) {
      double llhr = h1->GetBinCenter(ib+1);
      if ((llhr >= fMinLLHR) and (llhr <= fMaxLLHR)) {
        h1->SetBinContent(ib+1,1);
        h1->SetBinError  (ib+1,0);
      }
    }
    double total = h1->Integral();
    h1->Scale(1./total);
  }
  else {
//-----------------------------------------------------------------------------
// general case: kinematic info is available
// for each NObs need to generate multiple distributions - B(k)+S(NObs-k)
// generate momentum, background hypothesis
//-----------------------------------------------------------------------------
    TObjArray* as  = fHist.fLogLhs [NObs];
    TObjArray* ab  = fHist.fLogLhb [NObs];
    TObjArray* ar  = fHist.fLogLhr [NObs];
    TObjArray* arR = fHist.fLogLhrR[NObs];
//-----------------------------------------------------------------------------
// loop over the background configurations
//-----------------------------------------------------------------------------
    for (int nb=0; nb<=NObs; nb++) {
      TH1D* hs  = (TH1D*) as->At(nb);
      TH1D* hb  = (TH1D*) ab->At(nb);
      TH1D* hr  = (TH1D*) ar->At(nb);
      TH1D* hrR = (TH1D*) arR->At(nb);
//-----------------------------------------------------------------------------
// pseudoexperiments
//-----------------------------------------------------------------------------
      int    ns = NObs-nb;
      for (int ipe=0; ipe<NPe; ipe++) {
        double tot_lhb =  1;
        double tot_lhs =  1; 
        double p       = -1;
                                        // first: nb background events
        for (int i=0; i<nb; i++) {
          p         = bgr_mom();
          
          double lhb  = lh_bgr(p);           // 
          tot_lhb    *=lhb;

          double lhs = lh_sig(p);
          tot_lhs   *= lhs;
        }
                                        // next: ns=NObs-nb signal events
        for (int i=0; i<ns; i++) {
          p           = sig_mom();
          
          double lhb  = lh_bgr(p);           // 
          tot_lhb    *=lhb;

          double lhs = lh_sig(p);
          tot_lhs   *= lhs;
        }
          
        double llhs = log(tot_lhs);
        double llhb = log(tot_lhb);
    
        double llhr = llhb-llhs;

        if (fDebug.fRun != 0) {
          printf("tot_lhs, tot_lhb, llhs, llhb, llhr: %12.5e %12.5e %12.5e  %12.5e %12.5e \n",
                 tot_lhs, tot_lhb, llhs, llhb, llhr);
        }

        hs->Fill(llhs);
        hb->Fill(llhb);
        hr->Fill(llhr);
                                        // reduced likelihood - [0,1]
        
        double llhrR = ((llhr/NObs)-fMinLLHR)/(fMaxLLHR-fMinLLHR);
    
        hrR->Fill(llhrR);
        // if (fDebug.fRun > 0) {
        //                                 // fill histograms for the generated momentum
        //   fHist.gen_pbgr->Fill(p,1);
        //   fHist.gen_psig->Fill(p,1);
        // }
      }
      double total = hrR->Integral();
      hrR->Scale(1./total);
    }
  }
//-----------------------------------------------------------------------------
// fHist.log_lhrR is used to build the acceptance interval, normalize to one
//-----------------------------------------------------------------------------
  return 0;
}

//-----------------------------------------------------------------------------
// the total number of histograms is large, so organize them by directories
// N(
//-----------------------------------------------------------------------------
int TKinLH::save_hist(const char* Filename, const char* Option) {

  TFile* f = TFile::Open(Filename,Option);

  for (int nt=0; nt<MaxNx; nt++) {
    f->mkdir(Form("%05i",nt));
  }
                                        // current directory is still //
  for (int nobs=0; nobs<MaxNx; nobs++) {
    f->cd(Form("//%05i",nobs));

    TObjArray* as  = fHist.fLogLhs [nobs];
    TObjArray* ab  = fHist.fLogLhb [nobs];
    TObjArray* ar  = fHist.fLogLhr [nobs];
    TObjArray* arR = fHist.fLogLhrR[nobs];

    as ->Write(); //Form("llhs_%02i" ,nobs));
    ab ->Write(); //Form("llhb_%02i" ,nobs));
    ar ->Write(); //Form("llhr_%02i" ,nobs));
    arR->Write(); //Form("llhrR_%02i",nobs));
  }

  f->cd("//");

  fHist.prob_sig->Write();
  fHist.prob_bgr->Write();

  //  f->Write();
  f->Close();
  
  delete f;
  
  return 0;
}
  
//-----------------------------------------------------------------------------
// assume S in 
//-----------------------------------------------------------------------------
int TKinLH::test_coverage(double MuB, double SMin, double SMax, int NPoints) {
  return 0;
}
//-----------------------------------------------------------------------------
void TKinLH::Print(const char* Option) const {
  printf("%-20s: PMin: %10.3f PMax:%10.3f\n",GetName(), pmin, pmax);
}

}
