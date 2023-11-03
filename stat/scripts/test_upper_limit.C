///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "Stntuple/stat/upper_limit.hh"

using namespace stntuple ;

//-----------------------------------------------------------------------------
// tests
//-----------------------------------------------------------------------------
upper_limit* test_ul_001(double MuB, double MuS = 1, double CL = 0.9, int NObs = -1, const char* Name = "test_ul_001") {

  upper_limit* cg = new upper_limit(Name,CL);

  cg->construct_interval(MuB,MuS,NObs);
                                        // make histograms
  cg->make_prob_hist();

  char cname[100];
  sprintf(cname,"c_%s",Name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) c = new TCanvas(cname,Name,1200,800);
  else              c->cd();

  cg->fHist.fProb->Draw("text");
  cg->fHist.fInterval->Draw("sames");

  printf("IxMin, IxMax, Sum(P) : %3i %3i %12.5e\n",cg->fIxMin,cg->fIxMax,cg->fSump);

  return cg;
}

//-----------------------------------------------------------------------------
// construct belt
//-----------------------------------------------------------------------------
upper_limit* test_ul_002(double MuB = 0, double CL = 0.9, int NObs = -1, const char* Name = "test_ul_002") {

  upper_limit* cg = new upper_limit(Name,CL);

                                        // make histograms
  cg->construct_belt(MuB,0,35,35001,NObs);
  cg->make_belt_hist();

  char cname[100];
  sprintf(cname,"c_%s",Name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) c = new TCanvas(cname,Name,1200,800);
  else              c->cd();

  cg->belt(NObs)->fStack->Draw();

  return cg;
}

//-----------------------------------------------------------------------------
// test_cg_003: test 90% coverage
//-----------------------------------------------------------------------------
upper_limit* test_ul_003(double MuB, double SMin, double SMax, int NPoints, double CL = 0.9, const char* Name = "test_ul_003") {

  upper_limit* cg = new upper_limit(Name,CL);

  cg->SetNExp(1000);
  cg->fDebug.fTestCoverage = 1;
                                        // test itself
  
  cg->test_coverage(MuB,SMin,SMax,NPoints);
  
                                        // make histograms
  char cname[100];
  sprintf(cname,"c_%s",Name);
  
  TCanvas* c = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(cname);
  if (c == nullptr) {
    c = new TCanvas(cname,Name,1200,800);
  }

  cg->fHist.fCoverage->Draw("alp");
  
  return cg;
}

