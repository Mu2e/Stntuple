///////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
// 0  : all events
// 1  : passed events
// 2  : rejected events
// 
// 3  : events with N>=2 "good" particles, good = (p>80) and (Nhits>20)
// 4  : events with NHitsTF > 1
// 5  : UNUSED
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

#include "Stntuple/ana/TGenAnaModule.hh"

ClassImp(stntuple::TGenAnaModule)

namespace stntuple {
//-----------------------------------------------------------------------------
TGenAnaModule::TGenAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode       = 11;
  fGeneratorCode = 28;
}

//-----------------------------------------------------------------------------
TGenAnaModule::~TGenAnaModule() {
}


//-----------------------------------------------------------------------------
void TGenAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  GenpHist_t* hist = Hist;

  HBook1F(hist->fPdgCode,"pdg_code",Form("%s: PDG code"      ,Folder), 200,-1000, 1000,Folder);
  HBook1F(hist->fGenCode,"gen_code",Form("%s: generator code",Folder), 100, -10, 90,Folder);
  HBook1F(hist->fMomentum,"mom"    ,Form("%s: Momentum"      ,Folder), 200, 0, 200,Folder);
}


//-----------------------------------------------------------------------------
void TGenAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //     char name [200];
  //     char title[200];
  SimpHist_t* hist = Hist;

  HBook1F(hist->fPdgCode   ,"pdg_code",Form("%s: PDG code"     ,Folder), 200,-1000, 1000,Folder);
  HBook1F(hist->fMomentum  ,"mom"    ,Form("%s: Momentum"      ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fNStrawHits,"nsh"    ,Form("%s: N(straw hits"  ,Folder), 200, 0, 200,Folder);
  HBook1F(hist->fProperTime,"tprop"  ,Form("%s: proper time"   ,Folder), 200, 0,  50,Folder);
}


//-----------------------------------------------------------------------------
void TGenAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  HBook1F(hist->fEventNumber,"evtnum"  ,Form("%s: Event Number" ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum"  ,Form("%s: Run   Number" ,Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNGenp      ,"ngenp"   ,Form("%s: N(genp)"      ,Folder), 100 , 0,  100 ,Folder);
  HBook1F(hist->fNSimp[0]   ,"nsimp_0" ,Form("%s: N(simp)[all]" ,Folder), 100 , 0,  100 ,Folder);
  HBook1F(hist->fNSimp[1]   ,"nsimp_1" ,Form("%s: N(simp)[good]",Folder), 100 , 0,  100 ,Folder);
}

//_____________________________________________________________________________
void TGenAnaModule::BookHistograms() {

  //  char name [200];
  //  char title[200];

  TFolder* fol;
  TFolder* hist_folder;
  char     folder_name[200];

  DeleteHistograms();
  hist_folder = (TFolder*) GetFolder()->FindObject("Hist");
//-----------------------------------------------------------------------------
// book event histograms
//-----------------------------------------------------------------------------
  TString*  event_selection   [kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) { event_selection[i] = NULL; }

  event_selection[ 0] = new TString("all events");		// all events

  for (int i=0; i<kNEventHistSets; i++) {
    if (event_selection[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,event_selection[i]->Data());
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book GENP histograms
//-----------------------------------------------------------------------------
  TString*  genp_selection   [kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) { genp_selection[i] = NULL; }

  genp_selection[0] = new TString("all genparticles");

  for (int i=0; i<kNGenpHistSets; i++) {
    if (genp_selection[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,genp_selection[i]->Data());
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book SIMP histograms
//-----------------------------------------------------------------------------
  TString*  simp_selection   [kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) { simp_selection[i] = nullptr; }

  simp_selection[0] = new TString("all simparticles");
  simp_selection[1] = new TString("simparticles nsh>20");

  for (int i=0; i<kNSimpHistSets; i++) {
    if (simp_selection[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,simp_selection[i]->Data());
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
}

//-----------------------------------------------------------------------------
void TGenAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Part) {

  TLorentzVector p;

  Part->Momentum(p);

  GenpHist_t* hist = (GenpHist_t*) Hist;
  
  hist->fPdgCode->Fill(Part->GetPdgCode());
  hist->fGenCode->Fill(Part->GetStatusCode());
  hist->fMomentum->Fill(p.P());
}

//-----------------------------------------------------------------------------
void TGenAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Part) {

  const TLorentzVector* p;

  p = Part->StartMom();

  SimpHist_t* hist = (SimpHist_t*) Hist;
  
  hist->fPdgCode->Fill(Part->PDGCode());
  hist->fMomentum->Fill(p->P());
  hist->fNStrawHits->Fill(Part->NStrawHits());

  float t0 = Part->StartProperTime();
  float t1 = Part->EndProperTime  ();
  float dt = (t0 >= 0) ? t1-t0 : -1;

  hist->fProperTime->Fill(dt);
}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TGenAnaModule::FillEventHistograms(EventHist_t* Hist) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;

  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fEventNumber->Fill(event_number);
  hist->fRunNumber->Fill(run_number);
  hist->fNGenp->Fill(fNGenp);
  hist->fNSimp[0]->Fill(fNSimp[0]);
  hist->fNSimp[1]->Fill(fNSimp[1]);

}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TGenAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("GenpBlock" ,"TGenpBlock" ,&fGenpBlock);
  RegisterDataBlock("SimpBlock" ,"TSimpBlock" ,&fSimpBlock);
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}


//_____________________________________________________________________________
void TGenAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);
//-----------------------------------------------------------------------------
// GENP histograms
//-----------------------------------------------------------------------------
  TGenParticle* genp;

  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);
  }
//-----------------------------------------------------------------------------
// SIMP histograms
//-----------------------------------------------------------------------------
  TSimParticle* simp;

  for (int i=0; i<fNSimp[0]; i++) {
    simp = fSimpBlock->Particle(i);
    FillSimpHistograms(fHist.fSimp[0],simp);
    if (simp->NStrawHits() > 20) FillSimpHistograms(fHist.fSimp[1],simp);
  }
}



//_____________________________________________________________________________
int TGenAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
int TGenAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
// if there are several hits, use the first one
//-----------------------------------------------------------------------------
  fNGenp      = fGenpBlock->NParticles();
  fNSimp[0]   = fSimpBlock->NParticles();

  fNSimp[1]   = 0;
  for (int i=0; i<fNSimp[0]; i++) {
    TSimParticle* simp = fSimpBlock->Particle(i);
    const TLorentzVector* p = simp->StartMom();
    
    if ((p->P() > 80) && (simp->NStrawHits() > 20)) {
      fNSimp[1] += 1;
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TGenAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 3 : events with N(good particles) <= 2
//-----------------------------------------------------------------------------
  if ((GetDebugBit(3) == 1) and (fNSimp[1] >= 2)) {
    GetHeaderBlock()->Print(Form("N(simp[1]) = %i",fNSimp[1]));
  }
//-----------------------------------------------------------------------------
// bit 4: events with NHitsTF > 1
//-----------------------------------------------------------------------------
  if (GetDebugBit(4) == 1) {
    GetHeaderBlock()->Print(Form("NHits(TF) = %5i",fNGenp));
  }
}

//_____________________________________________________________________________
int TGenAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TGenAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

}
