///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

#include "Stntuple/stat/TKinLH.hh"

ClassImp(stntuple::TKinLH)

namespace stntuple {

//-----------------------------------------------------------------------------
TKinLH::TKinLH(const char* Name, double CL, double PMin, double PMax, int Debug) : TBelt(Name, CL) {

  fDebug.fRun           = 0;
  fDebug.fConstructBelt = 0;
  fDebug.fTestCoverage  = 0;
  fDebug.fMuMin         = 0;
  fDebug.fMuMax         = -1;
  pmin                  = PMin;
  pmax                  = PMax;
  fColor                = kBlue+2;
  
  fHist.fProb       = nullptr;
  fHist.fLlh        = nullptr;
  fHist.fInterval   = nullptr;
  fHist.fBelt       = nullptr;
  fHist.fBeltLo     = nullptr;
  fHist.fBeltHi     = nullptr;
  fHist.fBeltSp     = nullptr;
  fHist.fCoverage   = nullptr;

  init();
}

//-----------------------------------------------------------------------------
TKinLH::~TKinLH() {
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
  fHist.gen_pbgr = new TH1D(Form("h_%s_gen_pbgr",GetName()),title,200,100,110);
    
  title = Form("Generated P(sig) name:%s pmin=%5.1f",GetName(),pmin);
  fHist.gen_psig = new TH1D(Form("h_%s_gen_psig",GetName()),title,200,100,110);
  
  int     nbins_llh(20000), nbins_llhr(20000), nbins_llhrR(10000);
  
  double  llh_min  (-100), llh_max  (100);
  double  llhr_min (-100), llhr_max (100);
  double  llhrR_min(  -2), llhrR_max(  3);

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
      h      = new TH1D(Form("h_%s_llhs_%02i_%02i",GetName(),nb,ns),title.Data(),nbins_llh,llh_min,llh_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhs[ntot]->Add(h);

      title  = Form("Log(Lhb) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_%s_llhb_%02i_%02i",GetName(),nb,ns),title.Data(),nbins_llh,llh_min,llh_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhb[ntot]->Add(h);

      title  = Form("Log(Lhr) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_%s_llhr_%02i_%02i",GetName(),nb,ns),title.Data(),nbins_llhr,llhr_min,llhr_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhr[ntot]->Add(h);

      title  = Form("Log(LhrR) name:%s nb:%02i ns:%02i",GetName(),nb,ns);
      h      = new TH1D(Form("h_%s_llhrR_%02i_%02i",GetName(),nb,ns),title.Data(),nbins_llhrR,llhrR_min,llhrR_max);
      h->SetMarkerStyle(6);
      h->SetMarkerColor(fColor);
      h->SetLineColor  (fColor);
      fHist.fLogLhrR[ntot]->Add(h);
    }
  }

  for (int ix=0; ix<MaxNx; ix++) {
    name  = Form("h_%s_log_lhrR_1_n%02i",GetName(),ix);
    title = Form("%s_log_lhrR_1_n%02i",GetName(),ix);
    
    fHist.fLogLhrR_1[ix] = new TH1D(name,title,nbins_llhrR,llhrR_min,llhrR_max);

    name  = Form("h_%s_log_lhrR_2_n%02i",GetName(),ix);
    title = Form("%s_log_lhrR_2_n%02i",GetName(),ix);
    
    fHist.fLogLhrR_2[ix] = new TH1D(name,title,nbins_llhrR,-50,50);
  }

  name  = Form("h_%s_sum_log_lhrR_1",GetName());
  title = Form("%s_sum_log_lhrR_1",GetName());
  
  fHist.fSumLogLhrR_2 = new TH1D(name,title,nbins_llhrR,-50,50);
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
// three parameters, to maintain uniform interface
// Nobs is used only to determine the binomial probabilities
// the step can be completed only after all log_lhrR histograms are filled
// assume the histograms are read in
//-----------------------------------------------------------------------------
int TKinLH::construct_interval(double MuB, double MuS, int NObs) {

                                        // need a different function to initialize probability
                                        // coefficiencts for a given NObs
  double pb[MaxNx];

  int rc = init_truncated_poisson_dist(MuB,NObs,pb);
  if (rc < 0) return rc;
//-----------------------------------------------------------------------------
// next: for given MuB and MuS, construct LogLhrR_N histograms 
//-----------------------------------------------------------------------------
  double exp_mus = TMath::Exp(-MuS);
  TH1D*  h0      = (TH1D*) fHist.fLogLhrR[0]->At(0);
  int    nx      = h0->GetNbinsX();
  
  for (int nt=0; nt<MaxNx; nt++) {
    TObjArray* arR = fHist.fLogLhrR[nt];
    fHist.fLogLhrR_1[nt]->Reset();
    for (int nb=0; nb<=nt; nb++) {
      int    ns = nt-nb;
      double ps = exp_mus*pow(MuS,ns)/fFactorial[ns];
      
                                        // this is the absolute normalization of the corresponding histogram
      double wt = pb[nb]*ps;
      TH1D* h   = (TH1D*) arR->At(nb);
                                        // logLhrR_N is normalized to the Poisson probability P(MuB,MuS,NObs)
                                        // just summing over all hists with the same NObs
      fHist.fLogLhrR_1[nt]->Add(h,wt);
    }
  }
//-----------------------------------------------------------------------------
// next: construct 2-sided likelihood hist. nobs=0 is a special case, keep int in mind
// find max bin, assume it corresponds to ix=1
// start from finding a bin with the max probability density - it has to be done here,
// as the result depends on MuB and MuS
//-----------------------------------------------------------------------------
  double pmax  = -1;
  
  for (int nt=0; nt<MaxNx; nt++) {
                                        // for n=0 all bins are filled with a constant
                                        // such that the sum would equal to 1
    TH1D* h1 = fHist.fLogLhrR_1[nt];
    
    for (int ix=0; ix<nx; ix++) {
      double p = h1->GetBinContent(ix+1);
      if (p > pmax) {
        pmax  = p;
      }
    }
  }
  
  printf("TKinLH::construct_interval: MuB, MuS, NObs, pmax = %12.5e %12.5e %3i %12.5e\n",
         MuB,MuS,NObs,pmax); 
//-----------------------------------------------------------------------------
// now, create uniformly normalized distributions, 2-sided
//-----------------------------------------------------------------------------
  for (int nt=0; nt<MaxNx; nt++) {
    double p0 = fProb[nt];
//-----------------------------------------------------------------------------    
// h1 is supposed to be normalized to an integral (sum of contens of all bins) of 1
//-----------------------------------------------------------------------------
    TH1D*  h1 = fHist.fLogLhrR_1[nt];

    for (int ib=0; ib<nx; ib++) { 
      double p     = p0*(h1->GetBinContent(ib+1)/pmax);
      double log_p = -log(p);
      double wt    = p0*h1->GetBinContent(ib+1);
      if (nt < (MuB+MuS)) log_p = -log_p;
      fHist.fLogLhrR_2[nt]->Fill(log_p,wt);
    }
  }

  for (int nt=0; nt<MaxNx; nt++) {
    fHist.fSumLogLhrR_2->Add(fHist.fLogLhrR_2[nt]);
  }
//-----------------------------------------------------------------------------
// last step: define the interval in the likelihod_ratio space. remember - it is two-sided
// everything starts from bin nx/2+1
// assume nb is an even number
//-----------------------------------------------------------------------------
  int ipmax   = nx/2+1;
  TH1D* h_sum = fHist.fSumLogLhrR_2;

  double sump = h_sum->GetBinContent(ipmax);
  
                                        // bins symmetric wrt zero correspond to the same probability density
  for (int i=1; i<ipmax-1; i++) {
    double p1 = h_sum->GetBinContent(ipmax+i);
    double p2 = h_sum->GetBinContent(ipmax-i);
    sump = sump+p1+p2;
    if (sump >= fCL) {
                                        // done
      fInterval.fLlhrMin = 0;
      fInterval.fLlhrMax = h_sum->GetBinCenter(ipmax+i);  // interval bound - always positive
      fInterval.fProbTot = sump;
      fInterval.fPMax    = pmax;
      break;
    }
  }
  
  return 0;
}
  
//-----------------------------------------------------------------------------
// vary signal from SMin to SMax in NPoints, construct FC belt, fill belt histogram
// fBelt is the FC belt histogram
// avoid multiple useless re-initializations
//-----------------------------------------------------------------------------
int TKinLH::construct_belt(double MuB, double SMin, double SMax, int NPoints, int NObs) {

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
    double mus   = SMin+i*fBelt.fDy;

    int rc       = construct_interval(MuB,mus,NObs);
    // double lhmax = fInterval.fLlhrMax;
    // double lhmin = fInterval.fLlhrMin;

    if (rc == 0) {
      // double llh_lo = lhmin;
      // double llh_hi = lhmax;
      
      fBelt.fLlhInterval[5*i  ] = fInterval.fLlhrMin;
      fBelt.fLlhInterval[5*i+1] = fInterval.fLlhrMax;
      fBelt.fLlhInterval[5*i+2] = fInterval.fProbTot;
      fBelt.fLlhInterval[5*i+3] = fInterval.fPMax;
      fBelt.fLlhInterval[5*i+4] = -1;
    }
  }

  if (fDebug.fConstructBelt > 0) {
    for (int i=0; i<NPoints; i++) {
      printf("i,fBelt.fLlhInterval[5*i  ],fBelt.fLlhInterval[5*i+1]: %12.5f %12.5f %12.5f\n",
             fBelt.fLlhInterval[5*i  ],fBelt.fLlhInterval[5*i+1], fBelt.fLlhInterval[5*i+2]);
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
    delete fHist.fBeltSp;
  }
  
  fHist.fBeltLo   = new TH1D(Form("h_belt_lo_%s",GetName()),
                             Form("TBeltLH LO MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltHi   = new TH1D(Form("h_belt_hi_%s",GetName()),
                             Form("TBeltLH HI MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  fHist.fBeltSp = new TH1D(Form("h_belt_sp_%s",GetName()),
                             Form("TBeltLH NO MuB = %10.3f CL = %5.2f Nobs:%3i",fMuB,fCL,fNObs),
                             fBelt.fLlhNPoints,fBelt.fSMin,fBelt.fSMax);

  for (int ix=0; ix<fBelt.fLlhNPoints; ix++) {
    fHist.fBeltLo ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix  ]);
    fHist.fBeltHi ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+1]);
    fHist.fBeltSp ->SetBinContent(ix+1,fBelt.fLlhInterval[5*ix+2]);
  }

  fHist.fBeltLo->GetXaxis()->SetTitle("#mu_{S}");
  fHist.fBeltLo->GetYaxis()->SetTitle("LLH");
  
  // fHist.fBeltLo->SetFillStyle(fBelt.fFillStyle);
  // fHist.fBeltLo->SetFillColor(fBelt.fFillColor);
  fHist.fBeltLo->SetLineColor(fBelt.fFillColor);

  fHist.fBeltHi->SetFillStyle(fBelt.fFillStyle);
  fHist.fBeltHi->SetFillColor(fBelt.fFillColor);
  fHist.fBeltHi->SetLineColor(fBelt.fFillColor);

  fHist.fBeltSp->SetMarkerStyle(20);

  fHist.fBelt = new THStack(Form("hs_%s",GetName()),fHist.fBeltHi->GetTitle());
  fHist.fBelt->Add(fHist.fBeltLo);
  fHist.fBelt->Add(fHist.fBeltHi);
}

//-----------------------------------------------------------------------------
int TKinLH::read_hist(const char* Filename) {

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
      (*as) [nb] = (TH1D*) f->Get(Form("//%05i/h_%s_llhs_%02i_%02i" ,nobs,GetName(),nb,ns));
      (*ab) [nb] = (TH1D*) f->Get(Form("//%05i/h_%s_llhb_%02i_%02i" ,nobs,GetName(),nb,ns));
      (*ar) [nb] = (TH1D*) f->Get(Form("//%05i/h_%s_llhr_%02i_%02i" ,nobs,GetName(),nb,ns));
      (*arR)[nb] = (TH1D*) f->Get(Form("//%05i/h_%s_llhrR_%02i_%02i",nobs,GetName(),nb,ns));
    }

    // fHist.fLogLhrR_1[nobs]->Read(Form("h_llhrR1_%02i" ,nobs));
    // fHist.fLogLhrR_2[nobs]->Read(Form("h_llhrR2_%02i" ,nobs));
  }

  fHist.prob_sig = (TH1F*) f->Get(Form("h_%s_prob_sig",GetName()));
  fHist.prob_bgr = (TH1F*) f->Get(Form("h_%s_prob_bgr",GetName()));

  // fHist.fSumLogLhrR_2->Read(Form("h_sum_llhrR2_%02i" ,nobs));

  // f->Close();
  
  // delete f;
  
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
                                        // reduces likelihood
        double llhrR = llhr/NObs;
    
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

    // fHist.fLogLhrR_1[nobs]->Write(Form("h_llhrR1_%02i",nobs));
    // fHist.fLogLhrR_2[nobs]->Write(Form("h_llhrR2_%02i",nobs));
  }

  f->cd("//");

  fHist.prob_sig->Write();
  fHist.prob_bgr->Write();

  // fHist.fSumLogLhrR_N->Write(Form("h_sum_llhrR_%02i",nobs));

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
