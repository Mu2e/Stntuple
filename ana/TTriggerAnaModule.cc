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
// 3  : UNUSED
// 4  : events with NProtonStrawHits >= 20
// 5  : UNUSED
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/ana/TTriggerAnaModule.hh"

#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
//#include "Stntuple/obj/TStrawHit.hh"
#include "Stntuple/val/stntuple_val_functions.hh"

using stntuple::HistBase_t;
using stntuple::TrackPar_t;
//------------------------------------------------------------------------------
#include "Stntuple/ana/TTriggerAnaModule.hh"

ClassImp(stntuple::TTriggerAnaModule)

namespace stntuple {
//-----------------------------------------------------------------------------
TTriggerAnaModule::TTriggerAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPdgCode     = 11;
  fProcessCode = -1;
  fMinTrigMom  = 0.;
}

//-----------------------------------------------------------------------------
TTriggerAnaModule::~TTriggerAnaModule() {
}


//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTimeClusterHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TimeClusterHist_t* hist = (TimeClusterHist_t*) Hist;

  HBook1F(hist->fEnergy    ,"e"     ,Form("%s: cluster energy",Folder), 200,    0,  200,Folder);
  HBook1F(hist->fNHits     ,"nhits" ,Form("%s: nhits"         ,Folder), 200,    0,  200,Folder);
  HBook1F(hist->fR         ,"r"     ,Form("%s: radius"        ,Folder), 200,    0, 1000,Folder);
  HBook1F(hist->fTime      ,"time"  ,Form("%s: time"          ,Folder), 200,    0, 2000,Folder);
  
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookHelixHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  HelixHist_t* hist = (HelixHist_t*) Hist;

  HBook1F(hist->fNHits      ,"nHits"      ,Form("%s: nHits"         ,Folder), 200,    0,  200,Folder);
  HBook1F(hist->fNComboHits ,"nComboHits" ,Form("%s: nComboHits"    ,Folder), 100,    0,  100,Folder);
  HBook1F(hist->fLambda     ,"lambda"     ,Form("%s: lambda"        ,Folder),1500,  -1500, 1500,Folder);
  HBook1F(hist->fD0         ,"d0"         ,Form("%s: d0"            ,Folder),1000,-1000, 1000,Folder);
  HBook1F(hist->fNLoops     ,"nLoops"     ,Form("%s: nLoops"        ,Folder), 25,    0,   25,Folder);
  HBook1F(hist->fP          ,"p"          ,Form("%s: momentum"      ,Folder), 500,    0,  500,Folder);
  HBook1F(hist->fPt         ,"pt"         ,Form("%s: tMomentum"     ,Folder), 500,    0,  500,Folder);
  HBook1F(hist->fT0         ,"T0"         ,Form("%s: T0"            ,Folder), 1000,    0,  2000,Folder);
  HBook1F(hist->fChi2XYNDof  ,"Chi2XY"    ,Form("%s: Chi2XY/NDof"   ,Folder),  100,    0,    20,Folder);
  HBook1F(hist->fChi2PhiZNDof,"Chi2PhiZ"  ,Form("%s: Chi2PhiZ/NDof",Folder),   100,    0,    20,Folder);
  HBook1F(hist->fChi2TZNDof ,"Chi2TZ"     ,Form("%s: Chi2TZ/NDof"   ,Folder),  100,    0,    20,Folder);
  HBook1F(hist->fSimpPDG1   ,"PDG"        ,Form("%s: PDG"     ,Folder), 1000,    -500,  500,Folder);
  HBook1F(hist->fMCP        ,"pMC"        ,Form("%s: MC momentum"   ,Folder), 500,    0,  500,Folder);
  HBook1F(hist->fDp         ,"deltaP"     ,Form("%s: Delta P (Reconstructed - MC)"     ,Folder), 500,    -500,  500,Folder);
  HBook1F(hist->fDpt        ,"deltaPt"    ,Form("%s: Delta Pt (Reconstructed - MC)"     ,Folder), 500,    -500,  500,Folder);
  HBook1F(hist->fTZSlope    ,"tzSlope"    ,Form("%s: TZ slope"      ,Folder), 2000,    -1,  1,Folder);
  HBook1F(hist->fTZSlopeSig ,"tzSlopeSig" ,Form("%s: TZ slope significance",Folder), 200,    0,  20,Folder);
  HBook1F(hist->fHitRatio   ,"hitRatio"   ,Form("%s: Hit ratio"     ,Folder), 100,    0,  2,Folder);
  HBook1F(hist->fPropDir    ,"propDir"    ,Form("%s: Dir of propagation"   ,Folder), 3,    -1.5,  1.5,Folder);
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTrackHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TrackHist_t* hist = (TrackHist_t*) Hist;

  HBook1F(hist->fP      ,"p"    ,Form("%s: momentum" ,Folder), 1100,    0,  110,Folder);
  HBook1F(hist->fNActive,"nactv",Form("%s: N(Active)",Folder),  200,    0,  200,Folder);
  HBook1F(hist->fChi2Dof,"chi2d",Form("%s: Chi2/DOF" ,Folder),  200,    0,   40,Folder);
  HBook1F(hist->fT0     ,"t0"   ,Form("%s: T0"       ,Folder), 2000,    0, 2000,Folder);
  HBook1F(hist->fD0     ,"d0"   ,Form("%s: D0"       ,Folder),  200, -200,  200,Folder);
  HBook1F(hist->fZ0     ,"z0"   ,Form("%s: z0"       ,Folder),  200,-2000, 2000,Folder);
  HBook1F(hist->fTanDip ,"tdip" ,Form("%s: TanDip"   ,Folder),  250, -2.5,  2.5,Folder);
  HBook1F(hist->fAlgMask,"amask",Form("%s: AlgMask"  ,Folder),   10,    0,   10,Folder);
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookTriggerHistograms(HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  TriggerHist_t* hist = (TriggerHist_t*) Hist;
  float  npot[20] = {0, 9333333.333  , 12333333.333, 14666666.667 ,17000000.000 ,19333333.333 ,
                     21666666.667 ,24000000.000 ,26333333.333 ,29000000.000 ,32000000.000 ,35333333.333 ,39000000.000 ,
                     43333333.333 ,48333333.333 , 54333333.333 ,62666666.667 , 74666666.667 ,97666666.667, 1e8 };


  // {0,2333333.333,3333333.333,4000000.000,4666666.667,5333333.333,6000000.000,6666666.667,7333333.333,8000000.000,8666666.667,9333333.333,10000000.000,10666666.667,
  // 11333333.333,12000000.000,12666666.667,13333333.333,14000000.000,14666666.667,15333333.333,16000000.000,16666666.667,17666666.667,18666666.667,19666666.667,
  //                    20666666.667,21666666.667,22666666.667,23666666.667,25000000.000,26333333.333,27666666.667,29333333.333,31000000.000,33000000.000,35333333.333,38000000.000,41333333.333,45666666.667,52000000.000,63000000.000, 1e8};
  HBook1F(hist->fBits        ,"bits"    ,Form("%s: fired trigger bits" ,Folder), 300,    0,  300,Folder);
  HBook1F(hist->fNPassedPaths,"npassed" ,Form("%s: N(passed paths)"    ,Folder), 100,    0,  100,Folder);
  HBook1F(hist->fInstLum,"instLum",Form("%s: nPOT/#mu bunch"  ,Folder), 19, npot,Folder);
  
  TString trigNames[5] = {"tprDeMHighPStopTarg", "tprDePHighPStopTarg", "cprDeMHighPStopTarg", "cprDePHighPStopTarg", "allHighPStopTarg"}; 

  for (int i=0; i<5; i++){
    HBook1F(hist->fNPassedPathInstLum[i],Form("nPathLumi%d",i),Form("%s: nPOT/#mu bunch %s"  ,Folder, trigNames[i].Data()), 19, npot,Folder);
  }
  
  HBook1F(hist->fNPassedAprInstLum[0],"AprPassedLumi0",Form("%s APR highP: nPOT/#mu bunch "  ,Folder), 19, npot,Folder);
  HBook1F(hist->fNPassedAprInstLum[1],"AprPassedLumi1",Form("%s APR lowP: nPOT/#mu bunch "  ,Folder), 19, npot,Folder);
  HBook1F(hist->fNPassedAprInstLum[2],"AprPassedLumi2",Form("%s APR helix: nPOT/#mu bunch "  ,Folder), 19, npot,Folder);
  HBook1F(hist->fNPassedCprInstLum,"AprCprPassedLumi",Form("%s: nPOT/#mu bunch "  ,Folder), 19, npot,Folder);
  HBook1F(hist->fNPassedTprInstLum,"AprTprPassedLumi",Form("%s: nPOT/#mu bunch "  ,Folder), 19, npot,Folder);

}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::BookEventHistograms(stntuple::HistBase_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];
  EventHist_t* hist = (EventHist_t*) Hist;

  float  npot[20] = {0, 9333333.333  , 12333333.333, 14666666.667 ,17000000.000 ,19333333.333 ,
                     21666666.667 ,24000000.000 ,26333333.333 ,29000000.000 ,32000000.000 ,35333333.333 ,39000000.000 ,
                     43333333.333 ,48333333.333 , 54333333.333 ,62666666.667 , 74666666.667 ,97666666.667, 1e8 };

    //{0,2333333.333,3333333.333,4000000.000,4666666.667,5333333.333,6000000.000,6666666.667,7333333.333,8000000.000,8666666.667,9333333.333,10000000.000,10666666.667,
    //  11333333.333,12000000.000,12666666.667,13333333.333,14000000.000,14666666.667,15333333.333,16000000.000,16666666.667,17666666.667,18666666.667,19666666.667,
    //                     20666666.667,21666666.667,22666666.667,23666666.667,25000000.000,26333333.333,27666666.667,29333333.333,31000000.000,33000000.000,35333333.333,38000000.000,41333333.333,45666666.667,52000000.000,63000000.000, 1e8};;
  HBook1F(hist->fInstLum[0],"instLum0",Form("%s: nPOT/#mu bunch"  ,Folder), 300, 0,  1.e8,Folder);
  HBook1F(hist->fInstLum[1],"instLum1",Form("%s: nPOT/#mu bunch"  ,Folder), 19, npot,Folder);
  HBook1F(hist->fEventNumber,"evtnum",Form("%s: Event Number"  ,Folder), 1000, 0,  1.e4,Folder);
  HBook1F(hist->fRunNumber  ,"runnum",Form("%s: Run   Number"  ,Folder), 1000, 0,  1.e6,Folder);
  HBook1F(hist->fNTimeClusters ,"ntcl",Form("%s: N(time clusters)",Folder),   20, 0,    20,Folder);
  HBook1F(hist->fNHelices ,"nhel",Form("%s: N(helices)",Folder),   10, 0,    10,Folder);
  HBook1F(hist->fNTracks    ,"ntrk"  ,Form("%s: N(tracks)"     ,Folder),   10, 0,    10,Folder);
  HBook1F(hist->fMcMom     ,"mc_mom"   ,Form("%s: MC Particle Momentum"            ,Folder),1000,  0,200,Folder);
  HBook1F(hist->fMcCosTh   ,"mc_costh" ,Form("%s: MC Particle Cos(Theta) Lab"      ,Folder),100,-1,1,Folder);
  HBook1F(hist->fDtHelix        ,"hDeltaT0"       ,Form("%s: Helix Delta T0"                         ,Folder),500,-500,500,Folder);
  HBook1F(hist->fDpHelix        ,"hDeltaP"       ,Form("%s: Helix Delta P"                         ,Folder), 500, -500, 500, Folder);
  HBook1F(hist->fDptHelix       ,"hDeltaPt"      ,Form("%s: Helix Delta Pt"                        ,Folder), 500, -500, 500, Folder);
  HBook1F(hist->fDradiusHelix   ,"hDeltaRadius"  ,Form("%s: Helix Delta Radius"                    ,Folder), 500, -500, 500, Folder);
  HBook1F(hist->fDd0Helix       ,"hDeltaD0"      ,Form("%s: Helix Delta abs(D0)"                    ,Folder),1000,-1000,1000, Folder);
  HBook1F(hist->fDrcenterHelix  ,"hDeltaRcenter"   ,Form("%s: Helix Delta Rcenter"                 ,Folder),1000,-1000,1000, Folder);
  HBook1F(hist->fDtTrack        ,"tDeltaT0"       ,Form("%s: Track Delta T0"                         ,Folder),500,-500,500,Folder);
  HBook1F(hist->fDpTrack        ,"tDeltaP"       ,Form("%s: Track Delta P"                         ,Folder),500,-500,500,Folder);
  HBook1F(hist->fDptTrack        ,"tDeltaPt"       ,Form("%s: Track Delta Pt"                         ,Folder),500,-500,500,Folder);
  HBook1F(hist->fDz0Track        ,"tDeltaZ0"       ,Form("%s: Track Delta Z0"                         ,Folder),300,-3000,3000,Folder);
  HBook1F(hist->fDd0Track        ,"tDeltaD0"       ,Form("%s: Track Delta abs(D0)"                         ,Folder),1000,-1000,1000,Folder);
  HBook1F(hist->fHelixPEff        ,"hPEff"       ,Form("%s: Helix P Efficiency"                         ,Folder),10,0,10,Folder);
  //HBook1F(hist->fHelixDeltaTEff   ,"hDeltaTEff"       ,Form("%s: Helix Delta T0 Efficiency"                         ,Folder),10,0,10,Folder);
  //HBook1F(hist->fHelixPEff        ,"hD0Eff"       ,Form("%s: Helix Delta D0 Efficiency"                         ,Folder),10,0,10,Folder);
}

//_____________________________________________________________________________
void TTriggerAnaModule::BookHistograms() {

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
  int book_event_histset[kNEventHistSets];
  for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

  book_event_histset[ 0] = 1;		// all events

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] != 0) {
      sprintf(folder_name,"evt_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fEvent[i] = new EventHist_t;
      BookEventHistograms(fHist.fEvent[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book time cluster histograms
//-----------------------------------------------------------------------------
  int book_tpeak_histset[kNTimeClusterHistSets];
  for (int i=0; i<kNTimeClusterHistSets; i++) book_tpeak_histset[i] = 0;

  book_tpeak_histset[ 0] = 1;		// all time clusters
  book_tpeak_histset[ 1] = 1;		// nhits>20

  for (int i=0; i<kNTimeClusterHistSets; i++) {
    if (book_tpeak_histset[i] != 0) {
      sprintf(folder_name,"tpeak_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTimeCluster[i] = new TimeClusterHist_t;
      BookTimeClusterHistograms(fHist.fTimeCluster[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book helix histograms
//-----------------------------------------------------------------------------
  int book_helix_histset[kNHelixHistSets];
  for (int i=0; i<kNHelixHistSets; i++) book_helix_histset[i] = 0;

  book_helix_histset[ 0] = 1;		// all helices
  book_helix_histset[ 1] = 1;           //helix with a twin helix within deltaT plus/minus 10 ns
  book_helix_histset[ 2] = 1;           //helix with a twin helix with deltaT > 40 ns
  book_helix_histset[ 3] = 1;           //helix with a twin helix with deltaT > 100ns
  book_helix_histset[ 4] = 1;           //helix with bad Chi2XY
  book_helix_histset[ 5] = 1;           //helix with bad Chi2PhiZ
  book_helix_histset[ 6] = 1;           //helix with p>=80 MeV/c

  book_helix_histset[ 10] = 1;		// all Apr helices
  book_helix_histset[ 11] = 1;		// all Apr helices
  book_helix_histset[ 12] = 1;		// all Apr helices
  book_helix_histset[ 13] = 1;		// all Apr helices


  for (int i=0; i<kNHelixHistSets; i++) {
    if (book_helix_histset[i] != 0) {
      sprintf(folder_name,"helix_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fHelix[i] = new HelixHist_t;
      BookHelixHistograms(fHist.fHelix[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all    events
  book_track_histset[  1] = 1;		// tracks |D0| < 200

  for (int i=0; i<kNTrackHistSets; i++) {
    if (book_track_histset[i] != 0) {
      sprintf(folder_name,"trk_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrack[i] = new TrackHist_t;
      BookTrackHistograms(fHist.fTrack[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book trigger histograms
//-----------------------------------------------------------------------------
  int book_trigger_histset[kNTriggerHistSets];
  for (int i=0; i<kNTriggerHistSets; i++) book_trigger_histset[i] = 0;

  book_trigger_histset[  0] = 1;		// all    events

  for (int i=0; i<kNTriggerHistSets; i++) {
    if (book_trigger_histset[i] != 0) {
      sprintf(folder_name,"trig_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fTrigger[i] = new TriggerHist_t;
      BookTriggerHistograms(fHist.fTrigger[i],Form("Hist/%s",folder_name));
    }
  }

}

//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillEventHistograms(HistBase_t* Hist, double Weight) {
//   double            cos_th, xv, yv, rv, zv, p;
//   TLorentzVector    mom;

  EventHist_t* hist = (EventHist_t*) Hist;
  int event_number = GetHeaderBlock()->EventNumber();
  int run_number   = GetHeaderBlock()->RunNumber();

  hist->fInstLum[0]->Fill(GetHeaderBlock()->InstLum());
  hist->fInstLum[1]->Fill(GetHeaderBlock()->InstLum());
  hist->fEventNumber->Fill(event_number,Weight);
  hist->fRunNumber->Fill(run_number,Weight);
  hist->fNTimeClusters->Fill(fNTimeClusters,Weight);
  hist->fNHelices->Fill(fNHelices[0]+fNHelices[1],Weight);
  hist->fNTracks->Fill(fNTracks,Weight);
  hist->fMcMom->Fill(fMcMom,Weight);
  hist->fMcCosTh->Fill(fMcCosTh,Weight);

  // Loop over number of helices of multihelix events
  for (int i=0; i<fNHelices[1]-1; i++) {
    TStnHelix* helix = fAprHelixBlock->Helix(i);
    
    for (int j = i+1; j<fNHelices[1]; j++){
      TStnHelix* helixN = fAprHelixBlock->Helix(j);
      
      // Calculate differences between each pair of helices
      float hDeltaT = std::abs(helixN->fT0 - helix->fT0);
      float deltaP = helixN->P() - helix->P();
      float deltaPt = helixN->Pt() - helix->Pt();
      float deltaRadius = helixN->fRadius - helix->fRadius;
      float deltaD0 = std::abs(helixN->fD0) - std::abs(helix->fD0);
      float deltaRCenter = helixN->fRCent - helix->fRCent;
      
      // Fill event histograms
      hist->fDtHelix->Fill(hDeltaT,Weight);
      hist->fDpHelix->Fill(deltaP,Weight);
      hist->fDptHelix->Fill(deltaPt,Weight);
      hist->fDradiusHelix->Fill(deltaRadius,Weight);
      hist->fDd0Helix->Fill(deltaD0,Weight);
      hist->fDrcenterHelix->Fill(deltaRCenter,Weight);

      // Saves index of helices with deltaT0 less than 10 ns
      if (hDeltaT < 10.){
	fHelixIndexGood.push_back(i);
	fHelixIndexGood.push_back(j);
      }

      // Saves index of helices with deltaT0 less than 100 ns 
      else if (hDeltaT < 100.){
	fHelixIndexBad.push_back(i);
	fHelixIndexBad.push_back(j);
      }

      // Saves index of helices with deltaT0 greater than 100 ns
      else {
	fHelixIndexTail.push_back(i);
	fHelixIndexTail.push_back(j);
      }

      // Flags on DeltaT0
      /* if (GetDebugBit(4) == 1) {
	if (deltaT >= 500) {
	  GetHeaderBlock()->Print(Form("deltaT = %10.3f", deltaT));
	}
      }
      if (GetDebugBit(5) == 1) {
	if (deltaT < 10) {
	  GetHeaderBlock()->Print(Form("deltaT = %10.3f", deltaT));
	}
      }
       if (GetDebugBit(6) == 1) {
	 if ((deltaT > 10) && (deltaT < 100)) {
	  GetHeaderBlock()->Print(Form("deltaT = %10.3f", deltaT));
	}
	}*/
    }
  }

  // Fill helix moment efficiency plot
  int passed = fTriggerBlock->PathPassed(1);
  if (passed){//this requires the events to have 2 MC particles in the tracker with N-hits>15 and P>60 MeV/c
    std::array<int,4> helixMomCut={80, 85, 90, 95};
    //std::array<int,3> helixDeltaTCut={100, 90, 80};
    //std::array<int,3> helixD0Cut={600, 500, 400};
    for (int i = 0; i<fNHelices[1]; i++){
      TStnHelix* helix = fAprHelixBlock->Helix(i);
      hist->fHelixPEff->Fill(0);
      for (int j = 0; j < 4; j++){
        if (helix->P() > helixMomCut[j]) hist->fHelixPEff->Fill(j+1);
      } 
    }
  }
  
  // Fill helix DeltaT and DeltaD0 plots
  /*for (int i = 0; i<fNHelices-1; i++){
    TStnHelix* helix = fHelixBlock->Helix(i);
    for (int j = i+1; j<fNHelices; j++){
      TStnHelix* helixN = fHelixBlock->Helix(j);
      float DeltaT = std::abs(helix->fT0 - helixN->fT0);
      float DeltaD0 = std::abs(helix->fD0) - std::abs(helixN->fD0);
      hist->fHelixDeltaTEff->Fill(0);
      hist->fHelixD0Eff->Fill(0);
      for (int k = 0; k<3; k++){
	if (DeltaT < helixDeltaTCut[k]) hist->fHelixDeltaTEff->Fill(k+1);
	if (DeltaD0 < helixD0Cut[k]) hist->fHelixD0Eff->Fill(k+1);
      }
    }
    }*/
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillHelixHistograms(HistBase_t* Hist, TStnHelix* Helix) {
  HelixHist_t* hist = (HelixHist_t*) Hist;

  hist->fPropDir->Fill(Helix->fPropDir);
  hist->fNHits->Fill(Helix->fNHits);
  hist->fNComboHits->Fill(Helix->fNComboHits);
  hist->fLambda ->Fill(Helix->fLambda);
  hist->fD0->Fill(Helix->fD0);
  hist->fNLoops->Fill(std::abs(Helix->fNLoops));
  hist->fP->Fill(Helix->P());
  hist->fPt->Fill(Helix->Pt());
  hist->fT0->Fill(Helix->fT0);
  hist->fTZSlope->Fill(Helix->TZSlope());
  hist->fTZSlopeSig->Fill(Helix->TZSlopeSig());
  hist->fHitRatio->Fill(Helix->HitRatio());

  //Chi2s
  hist->fChi2XYNDof->Fill(Helix->fChi2XYNDof);
  hist->fChi2PhiZNDof->Fill(Helix->fChi2PhiZNDof);
  hist->fChi2TZNDof->Fill(Helix->fChi2TZNDof);

  //MC info
  hist->fSimpPDG1->Fill(Helix->fSimpPDG1);
  hist->fMCP->Fill(Helix->SimpMom1().P());
  float deltaP = Helix->P() - Helix->SimpMom1().P();
  hist->fDp->Fill(deltaP);
  float deltaPt = Helix->Pt() - Helix->SimpMom1().Pt();
  hist->fDpt->Fill(deltaPt);  
  
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTimeClusterHistograms(HistBase_t* Hist, TStnTimeCluster* TPeak) {
  TimeClusterHist_t* hist = (TimeClusterHist_t*) Hist;

  hist->fNHits->Fill(TPeak->fNHits);
  hist->fEnergy->Fill(TPeak->fClusterEnergy);
  hist->fTime->Fill(TPeak->fClusterTime);

  double x = TPeak->fClusterX;
  double y = TPeak->fClusterY;
  double r = sqrt(x*x+y*y);
  hist->fR->Fill(r);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTrackHistograms(HistBase_t* Hist, TStnTrack* Track) {

  TLorentzVector  mom;
  int             itrk;
  TrackPar_t*     tp;

  TrackHist_t* hist = (TrackHist_t*) Hist;
					// pointer to local track parameters
  itrk = Track->Number();
  tp   = fTrackPar+itrk;

  float na = Track->NActive();

  hist->fP->Fill (tp->fP);		// corrected momentum in the first point
  hist->fNActive->Fill(na);
  hist->fChi2Dof->Fill(Track->fChi2/(na-5.));
  hist->fT0->Fill(Track->fT0);
  hist->fD0->Fill(Track->fD0);
  hist->fZ0->Fill(Track->fZ0);
  hist->fTanDip->Fill(Track->fTanDip);
  hist->fAlgMask->Fill(Track->AlgMask());
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::FillTriggerHistograms(HistBase_t* Hist) {
  TriggerHist_t* hist = (TriggerHist_t*) Hist;

  int nb = fTriggerBlock->Paths()->GetNBits();

  for (int i=0; i<nb; i++) {
    int passed = fTriggerBlock->PathPassed(i);
    if (passed) hist->fBits->Fill(i);
  }

  int n_passed_paths = fTriggerBlock->NPassedPaths();
  hist->fNPassedPaths->Fill(n_passed_paths);

  hist->fInstLum->Fill(GetHeaderBlock()->InstLum());
  int  trgiBitsHitghPStopTarg[4] = {100, 101, 150, 151};
  bool passedHighPStopTarg(false);
  for (int i=0; i<4; ++i){
    if (fTriggerBlock->PathPassed(trgiBitsHitghPStopTarg[i])){
      hist->fNPassedPathInstLum[i]->Fill(GetHeaderBlock()->InstLum());
      passedHighPStopTarg = true;
    }
  }
  if (passedHighPStopTarg) hist->fNPassedPathInstLum[4]->Fill(GetHeaderBlock()->InstLum());
  
  if (fTriggerBlock->PathPassed(180)){
      hist->fNPassedAprInstLum[0]->Fill(GetHeaderBlock()->InstLum());
  }

  if (fTriggerBlock->PathPassed(190)){
      hist->fNPassedAprInstLum[1]->Fill(GetHeaderBlock()->InstLum());
  }

  if (fTriggerBlock->PathPassed(195)){
      hist->fNPassedAprInstLum[2]->Fill(GetHeaderBlock()->InstLum());
  }
  //  InstLum
}

//_____________________________________________________________________________
void TTriggerAnaModule::FillHistograms() {

//-----------------------------------------------------------------------------
// time peak histograms
//
// TPEAK_0: all track seeds
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTimeClusters; i++) {
    TStnTimeCluster* tpeak = fTimeClusterBlock->TimeCluster(i);
    FillTimeClusterHistograms(fHist.fTimeCluster[0],tpeak);
    if (tpeak->NHits() > 20) {
      FillTimeClusterHistograms(fHist.fTimeCluster[1],tpeak);
    }
  }
//-----------------------------------------------------------------------------
// helix histograms
//-----------------------------------------------------------------------------
  int  trgiBitsHitghPStopTarg[4] = {100, 101, 150, 151};
  for (int i=0; i<fNHelices[0]; i++) {
    TStnHelix* helix = fHelixBlock->Helix(i);
    FillHelixHistograms(fHist.fHelix[0],helix);
    if (helix->fChi2XYNDof > 5.)
      FillHelixHistograms(fHist.fHelix[4],helix);
    if (helix->fChi2PhiZNDof > 5.)
      FillHelixHistograms(fHist.fHelix[5],helix);
    if (helix->P()>=80.)
      FillHelixHistograms(fHist.fHelix[6],helix);        

    //bool passedHighPStopTarg(false);
    for (int i=0; i<4; ++i){
      if (fTriggerBlock->PathPassed(trgiBitsHitghPStopTarg[i])){
        FillHelixHistograms(fHist.fHelix[6],helix);        
      }
    }
    // if (fTriggerBlock->PathPassed(100)) FillTrackSeedHistograms(fHist.fTrackSeed[30],tseed,fWeight);
    // if (fTriggerBlock->PathPassed(101)) FillTrackSeedHistograms(fHist.fTrackSeed[31],tseed,fWeight);
  }

  // for (int i=0; i<fNTrackSeeds[1]; i++) {
  //   TStnTrackSeed* tseed = fAprTrackSeedBlock->TrackSeed(i);
  //   FillTrackSeedHistograms(fHist.fTrackSeed[20],tseed,fWeight);
  //   if (fTriggerBlock->PathPassed(180)) FillTrackSeedHistograms(fHist.fTrackSeed[21],tseed,fWeight);
  //   if (fTriggerBlock->PathPassed(190)) FillTrackSeedHistograms(fHist.fTrackSeed[21],tseed,fWeight);
  // }

//-----------------------------------------------------------------------------
// track histograms
//
// TRK_0: all tracks
// TRK_0: tracks    |D0| < 200 mm
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTracks; i++) {
    TStnTrack* trk = fTrackBlock->Track(i);
    FillTrackHistograms(fHist.fTrack[0],trk);
    if (fabs(trk->fD0) < 200.) {
      FillTrackHistograms(fHist.fTrack[1],trk);
    }
  }
//-----------------------------------------------------------------------------
// event histograms
//
// EVT_0: all events
//-----------------------------------------------------------------------------
//Clear the content of fHelixIndexGood/Bad/Tail
  fHelixIndexGood.clear();
  fHelixIndexBad.clear();
  fHelixIndexTail.clear();

  FillEventHistograms(fHist.fEvent[0]);

  for (size_t i=0; i<fHelixIndexGood.size(); i++) {
    TStnHelix* helix = fAprHelixBlock->Helix(fHelixIndexGood[i]);
    FillHelixHistograms(fHist.fHelix[1],helix);
  }

  for (size_t i=0; i<fHelixIndexBad.size(); i++) {
    TStnHelix* helix = fAprHelixBlock->Helix(fHelixIndexBad[i]);
    FillHelixHistograms(fHist.fHelix[2],helix);
  }

  for (size_t i=0; i<fHelixIndexTail.size(); i++) {
    TStnHelix* helix = fAprHelixBlock->Helix(fHelixIndexTail[i]);
    FillHelixHistograms(fHist.fHelix[3],helix);
  }
  
//-----------------------------------------------------------------------------
// trigger histograms
//
// TRIG_0: all events
//-----------------------------------------------------------------------------
  FillTriggerHistograms(fHist.fTrigger[0]);
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TTriggerAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("GenpBlock"            , "TGenpBlock"          , &fGenpBlock       );
  RegisterDataBlock("TimeClusterBlockDeCpr", "TStnTimeClusterBlock", &fTimeClusterBlock);
  RegisterDataBlock("HelixBlockDe"         , "TStnHelixBlock"      , &fHelixBlock      );
  RegisterDataBlock("HelixBlockAprHighP"   , "TStnHelixBlock"      , &fAprHelixBlock      );
  RegisterDataBlock("TrackBlockAprTrg"     , "TStnTrackBlock"      , &fTrackBlock      );
  //RegisterDataBlock("ClusterBlock"         , "TStnClusterBlock"    , &fClusterBlock    );
  RegisterDataBlock("TriggerBlock"         , "TStnTriggerBlock"    , &fTriggerBlock    );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();

  return 0;
}

//_____________________________________________________________________________
int TTriggerAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}

//_____________________________________________________________________________
void TTriggerAnaModule::InitTrackPar() {
  for (int i=0; i<fNTracks; i++) {
    TStnTrack*  track = fTrackBlock->Track(i);
    TrackPar_t* tp    = fTrackPar+i;

     tp->fP = track->fP0;
  }
}

//_____________________________________________________________________________
int TTriggerAnaModule::Event(int ientry) {

  //  double                p;
  //  TLorentzVector        mom;

  fGenpBlock->GetEntry(ientry);
  fTimeClusterBlock->GetEntry(ientry);
  fHelixBlock->GetEntry(ientry);
  fAprHelixBlock->GetEntry(ientry);
  fTrackBlock->GetEntry(ientry);
  //fClusterBlock->GetEntry(ientry);
  fTriggerBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// MC generator info
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  int           pdg_code, generator_code;

  fParticle = NULL;
  fNGenp    = fGenpBlock->NParticles();

  TStntuple* stnt = TStntuple::Instance();

  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fProcessCode)) {
      fParticle = genp;
      break;
    }
  }

  fProcessCode = -1;
  fWeight      =  1;

  TLorentzVector    mom (1.,0.,0.,0);

  if (fParticle) {
    if ((fParticle->GetStatusCode() == 7) && (fParticle->GetPdgCode() == 11)) {
//-----------------------------------------------------------------------------
// flat electrons
//-----------------------------------------------------------------------------
      fProcessCode = 7;

      fParticle->Momentum(mom);
      double e     = mom.E();
      fWeight      = stnt->DioWeightAl(e);
    }
  }

  fMcMom       = mom.P();
  fMcCosTh     = mom.Pz()/fMcMom;
//-----------------------------------------------------------------------------
// 
//-----------------------------------------------------------------------------
  fNTimeClusters = fTimeClusterBlock->NTimeClusters(); // all
  fNHelices[0]   = fHelixBlock->NHelices(); // all
  fNHelices[1]   = fAprHelixBlock->NHelices(); // all
  fNTracks       = fTrackBlock->NTracks();

  InitTrackPar();

  fNGoodTracks = 0;
  for (int i=0; i<fNTracks; i++) {
    TStnTrack* trk = fTrackBlock->Track(i);
    if (fabs(trk->fD0) < 200.) {
      if (trk->fP > 100) {
	if ((trk->Chi2Dof() < 5) && (trk->NActive() >= 20)) {
	  fNGoodTracks += 1;
	}
      }
    }
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TTriggerAnaModule::Debug() {

//-----------------------------------------------------------------------------
// bit 4: events with NProtonStrawHits >= 20
//-----------------------------------------------------------------------------
  // if (GetDebugBit(4) == 1) {
  //   if (fNProtonStrawHits >= 20) {
  //     GetHeaderBlock()->Print(Form("NStrawHits = %5i, NProtonStrawHits = %5i",
  // 				   fNStrawHits,fNProtonStrawHits));
  //   }
  // }
}

//_____________________________________________________________________________
int TTriggerAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TTriggerAnaModule::Test001() {

  // mu2e::HexMap* hmap      = new mu2e::HexMap();

  // mu2e::HexLK hex_index(0,0);

  // for (int i=0; i<40; i++) {
  //   hex_index = hmap->lk(i);
  //   printf(" i,l,k = %5i %5i %5i\n",i,hex_index._l,hex_index._k);
  // }
}

}
