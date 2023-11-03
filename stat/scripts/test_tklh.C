///////////////////////////////////////////////////////////////////////////////
// tests of the TKinLH class
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/stat/TKinLH.hh"
using namespace stntuple;
//-----------------------------------------------------------------------------
// plot normalized reduced likelihood
//-----------------------------------------------------------------------------
TKinLH* test_tklh_001(const char* Name  = "test_mllh_001",
                      double      CL    = 0.9,
                      double      PMin  = 102,
                      double      PMax  = 105, 
                      int         NPExp = 1000000, 
                      int         NObs  = 1      ) {
 
  TString   title;

  TKinLH* x1 = new TKinLH(Name,CL,PMin,PMax);
  x1->run(NObs,NPExp);
  
  // make probability histogram

  TString hname = Form("h_%s_prob",Name);

  TH1D* h_prob = new TH1D(hname.Data(),Form("P(bgr) for N(obs):%i",NObs),20000,-20,0);

  // double tot =  x1->fHist.log_lhrR[NObs]->Integral();

  // int nx = ((TH1D*) x1->fHist.fLogLhrR[NObs]->At(0))->GetNbinsX();

  // for (int nb=0; nb<=NObs; nb++) {
  //   for (int i=0; i<nx; i++) {
  //     double p    = x1->fHist.log_lhrR[NObs]->GetBinContent(i+1)/tot;
  //     double logp = log(p);
    
  //     h_prob->Fill(logp,p);
  //   }
  // }

  return x1;
}

//-----------------------------------------------------------------------------
// generate distributions in LLHR for a range of Nobs and store them
//-----------------------------------------------------------------------------
TKinLH* test_tklh_002(const char* Name = "t2",
                      double      CL   = 0.9,  
                      double      PMin = 102.0,
                      double      PMax = 105.0, 
                      long int    NExp = 100000,
                      int         NMax = TKinLH::MaxNx) {
  TStopwatch t;
  TKinLH* x1 = new TKinLH(Name,CL,PMin,PMax);

  t.Start();
  
  for (int nobs=0; nobs<NMax; nobs++) {
    x1->run(nobs,NExp);
    t.Stop();
    printf("nobs = %2i real_time:%10.3f cpu_time: %10.3f\n",nobs,t.RealTime(),t.CpuTime());
    // t.Print();
    t.Continue();
  }
  t.Stop();
//-----------------------------------------------------------------------------
// save all histograms once in the end
//-----------------------------------------------------------------------------
  TString fn = Form("tklh_%s_%0.f_%.0f.root",x1->GetName(),PMin*10,PMax*10);
  
  TFile* f;

  x1->save_hist(fn.Data(),"recreate");

  return x1;
}


TH1D* h_bgr_llhr;
TH1D* h_sig_llhr;
//-----------------------------------------------------------------------------
// read histogram file and construct a confidence interval
// kinematic histograms do not know anything about CL, but do know about the object name...
//  'Name' defines the input file name 
//-----------------------------------------------------------------------------
TKinLH* test_tklh_003(const char* Name    ,
                      double      CL          ,
                      double      PMin = 102.0,
                      double      PMax = 105.0,
                      double      MuB  =   0.1,
                      double      MuS  =   4.5) {

  TKinLH* x1 = new TKinLH(Name,CL,PMin,PMax);

  TString fn = Form("/projects/mu2e/hist/stntuple_stat/tklh_%s_%.0f_%.0f.root",Name,PMin*10,PMax*10);

  x1->read_hist(fn.Data());

  //  x1->fDebug.fConstructInterval = 1;
  
  x1->construct_interval(MuB,MuS);

  x1->plot_interval();

  return x1;
}

//-----------------------------------------------------------------------------
// read histogram file, construct a confidence belt, plot it
// also plot a data likelihood (not sure what it is for NObs=0)
//-----------------------------------------------------------------------------
TKinLH* test_tklh_004(const char* Name,
                      double      CL      = 0.9,
                      double      PMin    = 102,
                      double      PMax    = 105,
                      double      MuB     = 0.1,
                      double      SMin    = 1,
                      double      SMax    = 5,
                      int         NPoints = 10,
                      int         NObs    = 0)
{
  double p[] = {103., 104., 103, 104, 103., 104};
  
  TKinLH* x1 = new TKinLH(Name,CL,PMin,PMax);

  TString fn = Form("/projects/mu2e/hist/stntuple_stat/tklh_%s_%.0f_%.0f.root",Name,PMin*10,PMax*10);

  x1->read_hist(fn.Data());

  x1->construct_belt(MuB,SMin,SMax,NPoints,NObs,p);
  x1->make_belt_hist();

  char cname[100];
  sprintf(cname,"c_%s",Name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) c = new TCanvas(cname,Name,1200,800);
  else              c->cd();

  x1->fHist.fBeltHi->GetYaxis()->SetRangeUser(0,20);
  x1->fHist.fBelt->SetMaximum(15);
  x1->fHist.fBelt->Draw();
  x1->fHist.fBeltLhdt->Draw("same");
  // x1->fHist.fBeltNO1->Draw("same");
    
  return x1;
}


//-----------------------------------------------------------------------------
// assume belt is constructed, calculate LH
//-----------------------------------------------------------------------------
int test_tklh_005(TKinLH* X, double MuB=0.1, double MuS=1, int NObs = 0) {

  double p[] = {103.,104., 103, 104, 103., 104., 103., 104., 103., 104.};

  double llhrR;
  
  X->construct_interval(MuB,MuS);

  double wt = X->wt_data(MuB,MuS,NObs,p,&llhrR);

  double nlog_wt = -log(wt);
  
  TLine* line = new TLine(wt,1.e-2,wt,1.e-10);

  line->SetLineColor(kRed+2);
  line->Draw();

  printf("data: nobs, llhrR, wt, -log_wt : %i %10.3f %9.3e %9.3e\n",NObs,llhrR, wt, nlog_wt);
  return 0;
}


//-----------------------------------------------------------------------------
// test quick sort
//-----------------------------------------------------------------------------
void test_tklh_006(int NCalls = -1) {

  int npt = 6;
  TKinLH::sdata dt[] = {{2,-5.},{3,-5.},{1,2.},{0,4.},{2,-1.},{17,0}};

  TKinLH::fDebug_QuickSort = NCalls;
  
  TKinLH::quickSort(dt,0,npt-1);

  for (int i=0; i<npt; i++) {
    printf("i: %2i dt[i].bin: %2i dt[i].x: %10.3f\n",i,dt[i].bin,dt[i].x);
  }
}
