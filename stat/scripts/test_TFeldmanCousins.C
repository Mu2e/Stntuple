//-----------------------------------------------------------------------------
// it makes sens to have test script names including the name of the tested class explicitly
// example of use:
// --------------
// .L Stntuple/stat/scripts/test_TFeldmanCousins.C
// fc = new test_fc();
// fc->test_012(-1)   // plot median dis
//
// fc_test_001(double MuB, double MuS, double CL) : construct interval for a given (MuB,MuS)
// fc_test_002(double MuB, double CL)             : construct belt     for a given (MuS,MuB)
//
// fc_test_017(double MuB,0,10,1000,double CL)    : test 90% coverage
//-----------------------------------------------------------------------------

#include "TCanvas.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMatrixD.h"

#include "Stntuple/stat/TFeldmanCousins.hh"

using namespace stntuple;

class test_fc {
public:

  struct Hist_t {
    TGraph* fNEvents;
    TGraph* fSMean;
    TGraph* fSMedian;
    TGraph* fS90Med;
    TGraph* fExcl90Median;
    TH2D*   fBelt;
  } fHist;

  test_fc() {
    fHist.fNEvents = nullptr;
    fHist.fSMean   = nullptr;
    fHist.fSMedian = nullptr;
    fHist.fS90Med  = nullptr;
    fHist.fBelt    = nullptr;
    fHist.fExcl90Median = nullptr;
  }

  void test_009(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1);
					// experiment with the expected background  N
					// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
					// determine FC propability of a CL-level observation of a signal
					// use mean discovery
  void test_010(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1);
  
					// plot N(events) needed for discovery vs MuB in [0.,0.5] with step 0.001
  void test_011(double CL);
					// plot "median discovery" signal vs MuB in [0.,0.5] with step 0.001
  void test_012(double CL);
					// plot "mean discovery" signal vs MuB in [0.,0.5] with step 0.001
  void test_013(double CL);
					// plot "median 90%CL" excluded 
  void test_014();

					// plot 90% CL "exclusion" 
  void test_016(double CL);

};

//-----------------------------------------------------------------------------
// test_001: construct interval for a given (MuS,MuB)
//-----------------------------------------------------------------------------
stntuple::TFeldmanCousins* test_fc_001(double MuB, double MuS, double CL=0.9, int NObs = -1, const char* Name = "test_fc_001") {
  stntuple::TFeldmanCousins* fc(nullptr);

  if (fc == nullptr) fc = new stntuple::TFeldmanCousins(Name,CL);
  else               fc->SetCL(CL);

  fc->fDebug.fConstructInterval = 1;
  fc->ConstructInterval(MuB,MuS,NObs);
  fc->MakeProbHist();
  fc->PrintProbs(10);
  
  TH1D* h1 = (TH1D*) fc->fHist.fProb->Clone(Form("h1_interval_%s",fc->GetName()));;
  h1->Reset();
  h1->SetTitle("h1_interval");
  int nx = h1->GetNbinsX();
  
  for (int i=0; i<nx; i++) {
    if ((i>=fc->fIxMin) and (i<=fc->fIxMax)) h1->SetBinContent(i+1,fc->fHist.fProb->GetBinContent(i+1));
  }

  char cname[100];
  sprintf(cname,"c_%s",Name);
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) c = new TCanvas(cname,Name,1200,800);
  else              c->cd();

  double xmax = 4*(MuB+MuS);
  int i  = xmax / 5; 
  fc->fHist.fProb->GetXaxis()->SetRangeUser(0,i*5);
  fc->fHist.fProb->GetXaxis()->SetTitle("N(observed)");
  fc->fHist.fProb->SetTitle(Form("FC interval CL=%3.2f #mu_{B} = %4.3f #mu_S = %4.3f NObs=%2i",CL,MuB,MuS,NObs));
  fc->fHist.fProb->Draw("text+hist");

  h1->SetFillStyle(3005);
  h1->SetFillColor(kRed+2);
  h1->Draw("sames");

  return fc;
}

//-----------------------------------------------------------------------------
// fc_test_002: construct belt for a given (MuS,MuB)
//-----------------------------------------------------------------------------
stntuple::TFeldmanCousins* test_fc_002(double MuB, double CL, int NObs = -1, const char* Name = "test_fc_002") {
  stntuple::TFeldmanCousins* fc(nullptr);

  if (fc == nullptr) fc = new stntuple::TFeldmanCousins(Name,CL);
  else               fc->SetCL(CL);

  // fc->fDebug.fConstructBelt = 2;
  
  //  fc->ConstructBelt(MuB,0,35,3501/*001*/,NObs);
  fc->ConstructBelt(MuB,0,20,20001/*001*/,NObs);
  fc->MakeBeltHist();

  char cname[100];
  sprintf(cname,"c_%s",Name);

  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) c = new TCanvas(cname,Name,850,800);
  else              c->cd();

  fc->fHist.fBelt->SetTitle(Form("FC belt CL=%3.2f #mu_{B}=%4.3f NObs=%2i",CL,MuB,NObs));
  fc->fHist.fBelt->Draw();

  return fc;
}

//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine FC propability of a CL-level observation of a signal
//-----------------------------------------------------------------------------
void test_fc::test_009(double CL, double MuB, double  SMin, double SMax, int NPoints) {
  stntuple::TFeldmanCousins* fc(nullptr);

  if (fc == nullptr) fc = new stntuple::TFeldmanCousins("fc_009",CL);
  else               fc->SetCL(CL);

  fc->SetNExp(100000); // enough for defining 5

  double  x[NPoints], prob[NPoints];
  
  fc->DiscoveryProb(MuB,SMin,SMax,NPoints,x,prob);

  TCanvas* c = new TCanvas("c_fc_009","FC 009",1100,800);
  if (! c->GetShowEventStatus()) c->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_fc_009",int(MuB*1000)));
  gr->SetTitle(Form("FC discovery prob for bgr=%5.3f events",MuB));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");

  gPad->SetGridy(1);
  gPad->SetGridx(1);

  gPad->Modified();
  gPad->Update();
}

void test_fc::test_010(double CL, double MuB=0.1, double  SMin = 0, double SMax = 0, int NPoints = 1) {
//-----------------------------------------------------------------------------
// experiment with the expected background  N
// scan a range of signal hypotheses [SMin,SMax] with a step 'Step'
// determine FC propability of a CL-level observation of a signal
// use mean discovery
//-----------------------------------------------------------------------------
  stntuple::TFeldmanCousins* fc(nullptr);

  if (fc == nullptr) fc = new stntuple::TFeldmanCousins("fc010",CL);
  else               fc->SetCL(CL);

  fc->SetNExp(100000); // enough for defining 5

  double  x[NPoints], prob[NPoints];
  
  fc->DiscoveryProbMean(MuB,SMin,SMax,NPoints,x,prob);

  int i1 = -1;
  for (int i=0; i< NPoints; i++) {
    printf("i, x[i], prob[i] : %2i %10.3e %12.5e\n",i,x[i],prob[i]);
    if (prob[i] < 5.) continue;
    i1 = i-1;
    break;
  }

  printf("x1,y1, x2,y2, x3,y3: %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
  	 x[i1],prob[i1],x[i1+1],prob[i1+1],x[i1+2],prob[i1+2]);

  TMatrixD m(3,3);

  m(0,0) = x[i1  ]*x[i1  ]; m(0,1) = x[i1  ]; m(0,2) = 1;
  m(1,0) = x[i1+1]*x[i1+1]; m(1,1) = x[i1+1]; m(1,2) = 1;
  m(2,0) = x[i1+2]*x[i1+2]; m(2,1) = x[i1+2]; m(2,2) = 1;

  TMatrixD minv(TMatrixD::kInverted,m);

  TVectorD v1(3), v2(3);
  v1(0) = prob[i1  ];
  v1(1) = prob[i1+1];
  v1(2) = prob[i1+2];

  v2 = minv*v1;

  double a = v2(0);
  double b = v2(1);
  double c = v2(2);

  printf ("a,b,c = %12.5e %12.5e %12.5e\n",a,b,c);

  double xx = -b/(2*a) - sqrt(b*b/(4*a*a)-(c-5)/a);

  printf("xx = %12.5e\n",xx);
//-----------------------------------------------------------------------------
// have 3 points - draw a parabola, solve to 5
//-----------------------------------------------------------------------------
  TCanvas* cnv = new TCanvas("c_fc_010","FC 010",1100,800);
  if (! cnv->GetShowEventStatus()) cnv->ToggleEventStatus();
  gPad->SetCrosshair(1);

  TGraph* gr = new TGraph(NPoints,x,prob);

  gr->SetName(Form("gr_bgr_%06i_fc_010",int(MuB*1000)));
  gr->SetTitle(Form("FC discovery prob for bgr=%5.3f events",MuB));

  gr->SetMarkerStyle(20);
  gr->Draw("alp");

  gPad->SetGridy(1);
  gPad->SetGridx(1);

  gPad->Modified();
  gPad->Update();
}

//-----------------------------------------------------------------------------
// plot N(discovery) vs B
//-----------------------------------------------------------------------------
void test_fc::test_011(double CL) {

  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins("fc_011",CL);

  const char* name     = "c_fc_011";
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas("c_fc_PlotDiscoveryProbMean","c",1200,800);
  }

  int np(500);
  float x[np], y[np];
  
  for (int i=0; i<np; i++) {
    double mub = 0.001*(i+0.5);
    fc->ConstructInterval(mub,0);
    x[i] = mub;
    y[i] = fc->fIxMax;
  }

  if (fHist.fNEvents != nullptr) delete fHist.fNEvents;
  fHist.fNEvents = new TGraph(np,x,y);
  fHist.fNEvents->SetTitle("N(discovery events) vs #mu_{B}");
  fHist.fNEvents->Draw();
}

//-----------------------------------------------------------------------------
// plot "median discovery" signal vs B
//-----------------------------------------------------------------------------
void test_fc::test_012(double CL) {
  char name[100] = "fc_012";

  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins(name,CL);
  fc->SetNExp(100000);

  char cname[100];
  sprintf(cname,"c_%s",name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas(name,name,1200,800);
  }

  const int np(50); double x[np], y[np];

  double smin=1;
  double smax=11;
  
  const int ns  = 11; double s[ns], smean[ns];
  
  for (int i=0; i<np; i++) {
    double mub = 0.01*(i+0.5);
    fc->DiscoveryProb(mub,smin,smax,ns,s,smean);

    // solve for probability = 0.5

    x[i] = mub;
    fc->SolveFor(0.5,s,smean,ns,&y[i]);
  }

  if (fHist.fSMedian != nullptr) delete fHist.fSMedian;

  fHist.fSMedian = new TGraph(np,x,y);
  fHist.fSMedian->SetTitle("S(median) vs #mu_{B}");

  fHist.fSMedian->Draw();
}

//-----------------------------------------------------------------------------
// plot "mean discovery" signal vs B
//-----------------------------------------------------------------------------
void test_fc::test_013(double CL) {
  char name[100] = "fc_013";
  
  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins(name,CL);
  fc->SetNExp(100000);

  char cname[100];
  sprintf(cname,"c_%s",name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) {
    c = new TCanvas(cname,cname,1200,800);
  }

  const int np(50); double x[np], y[np];

  double smin=1;
  double smax=11;
  
  const int ns  = 11; double s[ns], smean[ns];
  
  for (int i=0; i<np; i++) {
    double mub = 0.01*(i+0.5);
    fc->DiscoveryProbMean(mub,smin,smax,ns,s,smean);

    // solve for 5

    x[i] = mub;
    fc->SolveFor(5.,s,smean,ns,&y[i]);
  }

  if (fHist.fSMean != nullptr) delete fHist.fSMean;

  fHist.fSMean = new TGraph(np,x,y);
  fHist.fSMean->SetName(Form("g_%s",name));;
  fHist.fSMean->SetTitle("S(mean) vs #mu_{B}");

  fHist.fSMean->Draw();
}

//-----------------------------------------------------------------------------
// plot 90% CL median vs vs B
//-----------------------------------------------------------------------------
void test_fc::test_014() {
  char name[100] = "fc_014";
  
  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins(name,0.9);
  fc->SetNExp(100000);

  char cname[100];
  sprintf(cname,"c_%s",name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas(name,name,1200,800);
  }

  const int np(100); double x[np], y[np];

  double smin=1;
  double smax=11;
  
  const int ns  = 1000;  double s[ns];
  const int nx  = fc->MaxNx;
  double  belt[nx][2];
  
  for (int i=0; i<np; i++) {
    double mub = 0.01*(i+0.5);
    // look for the level where belt[0][1] = 0
    //    printf(" -- i = %i mub = %10.5f\n",i,mub);
    fc->ConstructBelt(mub,smin,smax,ns);
    // find

    double yy = fc->fBelt.fSign[0][1];
    double step = (smax-smin)/(ns-1);

    //    printf(" --- yy, step = %12.5e %12.5e\n",yy,step);
    
//    fc->ConstructBelt(mub,yy-step/2,yy+step/2,ns,s,belt);

    //    printf(" --- limit: = %12.5e\n",belt[0][1]);

    double cl90 = fc->fBelt.fSign[0][1];

    x[i] = mub;
    y[i] = cl90;
  }

  if (fHist.fS90Med != nullptr) delete fHist.fS90Med;

  fHist.fS90Med = new TGraph(np,x,y);
  fHist.fS90Med->SetTitle("S90(median) vs #mu_{B}");

  fHist.fS90Med->Draw();
}


//-----------------------------------------------------------------------------
// plot 90% CL "exclusion" 
//-----------------------------------------------------------------------------
void test_fc::test_016(double CL) {
  char name[100] = "fc_test_016";

  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins(name,CL);
  fc->SetNExp(100000);

  char cname[100];
  sprintf(cname,"c_%s",name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(name);
  if (c == nullptr) {
    c = new TCanvas(name,name,1200,800);
  }

  const int np(20 /* 50 */ ); double x[np], y[np];

  double smin=1;
  double smax=11;
  
  const int ns  = 11; double s[ns], prob[ns];
  
  for (int i=0; i<np; i++) {
      //    double mub = 0.01*(i+0.5);
    double mub = 0.05*(i+0.5);

    double s1 = smin;
    double s2 = smax;

    while (s2-s1 > 0.01) {
      fc->UpperLimit(mub,s1,s2,ns,s,prob);

      int loc  = -1;
      for (int is=0; is<ns; is++) {
	if ((prob[is] < 0.5) and (prob[is+1] > 0.5)) {
	  loc = is;
	  break;
	}
      }

      s1 = s[loc];
      s2 = s[loc+1];
    }
    
    x[i] = mub;
    y[i] = (s1+s2)/2;
  }

  if (fHist.fExcl90Median != nullptr) delete fHist.fExcl90Median;

  fHist.fExcl90Median = new TGraph(np,x,y);
  fHist.fExcl90Median->SetTitle("Excl(median) vs #mu_{B}");
  fHist.fExcl90Median->Draw();
}


//-----------------------------------------------------------------------------
// fc_test_017: test 90% coverage
// fc17 = fc_test_017(0., 0., 10,1000,0.9,100000)
//-----------------------------------------------------------------------------
stntuple::TFeldmanCousins* test_fc_017(double MuB,
                                       double CL  ,
                                       double SMin,       
                                       double SMax,
                                       int    NPoints,
                                       int    NObs      = -1, // standard
                                       int    NPsExp    = -1,
                                       const char* Name = "test_fc_017",
                                       int DebugLevel   = -1) {

  stntuple::TFeldmanCousins* fc = new stntuple::TFeldmanCousins(Name,CL);

  if (NPsExp > 0) fc->SetNExp(NPsExp);
                                        // test itself
  if (DebugLevel > 0) fc->fDebug.fTestCoverage = DebugLevel;

  fc->TestCoverage(MuB,SMin,SMax,NPoints);

                                        // make histograms
  char cname[100];
  sprintf(cname,"c_%s",Name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(Name);

  if (c == nullptr) c = new TCanvas(cname,Name,1200,800);
  else              c->cd();

  fc->fHist.fCoverage->Draw("alp");
  
  return fc;
}
