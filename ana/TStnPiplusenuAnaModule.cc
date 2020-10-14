//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// Tmp(0) : nax seg
// Tmp(1) : nst seg
// 
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
// 
//  3  : events with set C tracks and 70mm < |dx|  < 90 mm
//  4  : events with tan(dip) > 1.4
///////////////////////////////////////////////////////////////////////////////
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TEnv.h"
#include "TSystem.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"
#include "Stntuple/alg/TStntuple.hh"
#include "Stntuple/geom/TDisk.hh"
#include "Stntuple/val/stntuple_val_functions.hh"
//------------------------------------------------------------------------------
// Mu2e offline includes
//-----------------------------------------------------------------------------
#include "Stntuple/ana/TStnPiplusenuAnaModule.hh"

ClassImp(TStnPiplusenuAnaModule)

//-----------------------------------------------------------------------------
TStnPiplusenuAnaModule::TStnPiplusenuAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  fPtMin  = 1.;
  fTrackNumber.Set(100);

  fFillDioHist     = 1;

  fMinT0 = 0; // do not cut on time by default

  fTrackID = new TStnTrackID();
  fLogLH   = new TEmuLogLH();
//-----------------------------------------------------------------------------
// MC truth: define which MC particle to consider as signal
//-----------------------------------------------------------------------------
  fPdgCode       = 11;
  fGeneratorCode = 2;			// conversionGun, 28:StoppedParticleReactionGun
}

//-----------------------------------------------------------------------------
TStnPiplusenuAnaModule::~TStnPiplusenuAnaModule() {
}


//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::BookGenpHistograms(GenpHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP      ,"p"       ,Form("%s: Momentum"     ,Folder),1000,     0, 200,Folder);
  HBook1F(Hist->fPdgCode[0],"pdg_code_0",Form("%s: PDG Code[0]"     ,Folder),200, -100, 100,Folder);
  HBook1F(Hist->fPdgCode[1],"pdg_code_1",Form("%s: PDG Code[1]"     ,Folder),500, -2500, 2500,Folder);
  HBook1F(Hist->fGenID  ,"gen_id"  ,Form("%s: Generator ID" ,Folder), 100,     0, 100,Folder);
  //HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  5400, 6400,Folder);
  HBook1F(Hist->fX0     ,"x0"      ,Form("%s: X0"           ,Folder), 50,  -4000, -3800,Folder);
  HBook1F(Hist->fY0     ,"y0"      ,Form("%s: Y0"           ,Folder), 50,  -150, 150,Folder);
  HBook1F(Hist->fZ0     ,"z0"      ,Form("%s: Z0"           ,Folder), 500,  3500, 5500,Folder);
  HBook1F(Hist->fT0     ,"t0"      ,Form("%s: T0"           ,Folder), 200,     0, 2000,Folder);
  HBook1F(Hist->fR0     ,"r"       ,Form("%s: R0"           ,Folder), 100,     0,  100,Folder);
  HBook1F(Hist->fCosTh  ,"cos_th"  ,Form("%s: Cos(Theta)"   ,Folder), 200,   -1.,   1.,Folder);
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::BookTrackHistograms(TrackHist_t* Hist, const char* Folder) {
//   char name [200];
//   char title[200];

  HBook1F(Hist->fP[0]       ,"p"        ,Form("%s: Track P(Z1)"       ,Folder), 400,  90  ,120. ,Folder);
  HBook1F(Hist->fP[1]       ,"p_1"      ,Form("%s: Track P(total)[1]" ,Folder), 100, 104.5,105.5,Folder);
  HBook1F(Hist->fP[2]       ,"p_2"      ,Form("%s: Track P(total)[1]" ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP0         ,"p0"       ,Form("%s: Track P(Z0)"       ,Folder),1000,   0  ,200. ,Folder);
  HBook1F(Hist->fP2         ,"p2"       ,Form("%s: Track P(z=-1540)"  ,Folder),1000,   0  ,200. ,Folder);
  HBook1D(Hist->fPDio       ,"pdio"     ,Form("%s: Track P(DIO WT)"   ,Folder), 400,  90  ,110. ,Folder);
  Hist->fPDio->Sumw2(kTRUE);

  HBook1F(Hist->fFitMomErr  ,"momerr"   ,Form("%s: Track FitMomError" ,Folder), 200,   0  ,  1. ,Folder);
  HBook1F(Hist->fPFront     ,"pf"       ,Form("%s: Track P(front)   " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFront    ,"dpf"      ,Form("%s: Track P-P(front) " ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront0   ,"dp0f"     ,Form("%s: Track P0-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fDpFront2   ,"dp2f"     ,Form("%s: Track P2-P(front)" ,Folder), 200,  -5. ,  5. ,Folder);
  HBook1F(Hist->fPStOut     ,"pstout"   ,Form("%s: Track P(ST_Out)  " ,Folder), 400,  90  ,110. ,Folder);
  HBook1F(Hist->fDpFSt      ,"dpfst"    ,Form("%s: Track Pf-Psto"     ,Folder), 200,  -5  ,  5. ,Folder);
  HBook2F(Hist->fDpFVsZ1    ,"dpf_vs_z1",Form("%s: Track DPF Vs Z1"   ,Folder), 200, -2000.,0,200,-5.,5,Folder);

  HBook1F(Hist->fPt         ,"pt"       ,Form("%s: Track Pt"          ,Folder), 600, 75,95,Folder);
  HBook1F(Hist->fCosTh      ,"costh"    ,Form("%s: Track cos(theta)"  ,Folder), 100,-1,1,Folder);
  HBook1F(Hist->fChi2       ,"chi2"     ,Form("%s: Track chi2 total"  ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fNDof       ,"ndof"     ,Form("%s: Number of DOF"     ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fChi2Dof    ,"chi2d"    ,Form("%s: track chi2/N(dof)" ,Folder), 500, 0, 10,Folder);
  HBook1F(Hist->fNActive    ,"nactv"    ,Form("%s: N(active)"         ,Folder), 200, 0,200,Folder);
  HBook1F(Hist->fT0         ,"t0"       ,Form("%s: track T0"          ,Folder), 200, 0,2000,Folder);
  HBook1F(Hist->fT0Err      ,"t0err"    ,Form("%s: track T0Err"       ,Folder), 100, 0,  10,Folder);
  HBook1F(Hist->fQ          ,"q"        ,Form("%s: track Q"           ,Folder),   4,-2,   2,Folder);
  HBook1F(Hist->fFitCons[0] ,"fcon"     ,Form("%s: track fit cons [0]",Folder), 200, 0,   1,Folder);
  HBook1F(Hist->fFitCons[1] ,"fcon1"    ,Form("%s: track fit cons [1]",Folder), 1000, 0,   0.1,Folder);
  HBook1F(Hist->fD0         ,"d0"       ,Form("%s: track D0      "    ,Folder), 200,-200, 200,Folder);
  HBook1F(Hist->fZ0         ,"z0"       ,Form("%s: track Z0      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fTanDip     ,"tdip"     ,Form("%s: track tan(dip)"    ,Folder), 200, 0.0 ,2.0,Folder);
  HBook1F(Hist->fResid      ,"resid"    ,Form("%s: hit residuals"     ,Folder), 500,-0.5 ,0.5,Folder);
  HBook1F(Hist->fAlgMask    ,"alg"      ,Form("%s: algorithm mask"    ,Folder),  10,  0, 10,Folder);

  HBook1F(Hist->fDt         ,"dt"       ,Form("%s: track delta(T)"    ,Folder), 200,-20  ,20 ,Folder);
  HBook1F(Hist->fChi2Match  ,"chi2tcm"  ,Form("%s: chi2(t-c match)"   ,Folder), 250,  0  ,250 ,Folder);

  HBook1F(Hist->fDx         ,"dx"       ,Form("%s: track delta(X)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDy         ,"dy"       ,Form("%s: track delta(Y)"    ,Folder), 200,-500 ,500,Folder);
  HBook1F(Hist->fDz         ,"dz"       ,Form("%s: track delta(Z)"    ,Folder), 200,-250 ,250,Folder);
  HBook1F(Hist->fDu         ,"du"       ,Form("%s: track-cluster DU)" ,Folder), 250,-250 ,250,Folder);
  HBook1F(Hist->fDv         ,"dv"       ,Form("%s: track-cluster DV)" ,Folder), 200,-100 ,100,Folder);
  HBook2F(Hist->fDvVsDu     ,"dv_vs_du" ,Form("%s: Track Dv Vs Du"    ,Folder), 100, -250,250,100,-100.,100,Folder);
  HBook1F(Hist->fPath       ,"path"     ,Form("%s: track sdisk"       ,Folder),  50,   0 ,500,Folder);
  HBook2F(Hist->fDuVsPath   ,"du_vs_path",Form("%s: Track Du Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDucVsPath  ,"duc_vs_path",Form("%s: T-C Duc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsPath   ,"dv_vs_path",Form("%s: T-C  Dv Vs Path"  ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDvcVsPath  ,"dvc_vs_path",Form("%s: T-C Dvc Vs Path" ,Folder),  50,   0 ,500,200,-200.,200.,Folder);
  HBook2F(Hist->fDtVsPath   ,"dt_vs_path",Form("%s: T-C DT Vs Path"   ,Folder),  50,   0 ,500,100,  -5.,  5.,Folder);
  HBook2F(Hist->fDuVsTDip   ,"du_vs_tdip",Form("%s: Track Du Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);
  HBook2F(Hist->fDvVsTDip   ,"dv_vs_tdip",Form("%s: Track Dv Vs TDip" ,Folder), 100, 0.5 ,1.5,200,-200.,200.,Folder);

  HBook1F(Hist->fZ1         ,"z1"       ,Form("%s: track Z1      "    ,Folder), 200,-2000,2000,Folder);
  HBook1F(Hist->fNClusters  ,"ncl"      ,Form("%s: track N(clusters)" ,Folder),  10, 0   , 10,Folder);
  HBook1F(Hist->fVaneID     ,"vid"      ,Form("%s: track vane ID"     ,Folder),  10,-5   ,  5,Folder);
  HBook1F(Hist->fXCal       ,"xcal"     ,Form("%s: track XCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYCal       ,"ycal"     ,Form("%s: track YCal"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZCal       ,"zcal"     ,Form("%s: track ZCal"        ,Folder), 200, 1500,3500,Folder);
  HBook1F(Hist->fXTrk       ,"xtrk"     ,Form("%s: track XTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fYTrk       ,"ytrk"     ,Form("%s: track YTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fRTrk       ,"rtrk"     ,Form("%s: track RTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fZTrk       ,"ztrk"     ,Form("%s: track ZTrk"        ,Folder), 200,-1000,1000,Folder);
  HBook1F(Hist->fECl        ,"ecl"      ,Form("%s: cluster E"         ,Folder), 300, 0   ,150,Folder);
  HBook1F(Hist->fEClEKin    ,"ecl_ekin" ,Form("%s: cluster E/Ekin(mu)",Folder), 200, 0   ,2,Folder);
  HBook1F(Hist->fEp         ,"ep"       ,Form("%s: track E/P"         ,Folder), 300, 0   ,1.5,Folder);
  HBook2F(Hist->fEpVsPath   ,"ep_vs_path",Form("%s: E/P Vs Path"      ,Folder),  50,   0 ,500,150,  0.,  1.5,Folder);
  HBook2F(Hist->fNHVsStation,"nh_vs_st" ,Form("%s: N(hits) Vs Station",Folder),  40, 0,40,10,-0.5,9.5,Folder);
  HBook2F(Hist->fNHVsNSt    ,"nh_vs_nst",Form("%s: N(hits) Vs NSt"    ,Folder),  10,-0.5,9.5,40,-0.5,39.5,Folder);

  HBook1F(Hist->fRSlope     ,"rslope"   ,Form("%s: Res Slope"         ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fXSlope     ,"xslope"   ,Form("%s: Res/Sig Slope"     ,Folder), 200,-20 , 20,Folder);

  HBook2F(Hist->fEpVsDt     ,"ep_vs_dt" ,Form("%s: E/P vs Dt"         ,Folder), 200, -10, 10,150,0.,1.5,Folder);
  HBook1F(Hist->fEleLogLHCal,"ele_llh_c",Form("%s: ELE Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fMuoLogLHCal,"muo_llh_c",Form("%s: MUO Log(LH) Cal"   ,Folder), 200,-100,  0,Folder);
  HBook1F(Hist->fLogLHRCal  ,"llhr_cal" ,Form("%s: LogLH(e/m) Cal"    ,Folder), 200,-100,100,Folder);
  HBook1F(Hist->fLogLHRDeDx ,"llhr_dedx",Form("%s: LogLH(e/m) De/Dx"  ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRXs   ,"llhr_xs"  ,Form("%s: LogLH(e/m) XSlope" ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHRTrk  ,"llhr_trk" ,Form("%s: LogLH(e/m) Trk"    ,Folder), 200,-20 , 20,Folder);
  HBook1F(Hist->fLogLHR     ,"llhr"     ,Form("%s: LogLH(e/m)"        ,Folder), 200,-100 ,100,Folder);

  HBook1F(Hist->fPdgCode    ,"pdg"      ,Form("%s: track PDG code"    ,Folder), 100,-50,50,Folder);
  HBook1F(Hist->fFrGH       ,"fgh"      ,Form("%s: Fraction Goog Hits",Folder), 100, 0,1,Folder);

  HBook2F(Hist->fNEPlVsNHPl ,"nep_vs_nhp",Form("%s: Track NEXP vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fNDPlVsNHPl ,"ndp_vs_nhp",Form("%s: Track NDIF vs NHit",Folder), 100, 0,100,100,0.,100,Folder);
  HBook2F(Hist->fChi2dVsNDPl,"chi2d_vs_ndp",Form("%s: Track Chi2/Dof vs NDP",Folder), 30, 0,30,100,0.,10,Folder);
  HBook2F(Hist->fDpFVsNDPl  ,"dpf_vs_ndp"  ,Form("%s: Track DpF vs NDP",Folder)     , 30, 0,30,100,-5,5,Folder);

  HBook1F(Hist->fFrE1   ,"fre1"   ,Form("%s: E1/Etot"       ,Folder),200, 0,  1,Folder);
  HBook1F(Hist->fFrE2   ,"fre2"   ,Form("%s: (E1+E2)/Etot"  ,Folder),200, 0,  1,Folder);
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::BookEventHistograms(EventHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fEleCosTh  ,"ce_costh" ,Form("%s: Conversion Electron Cos(Theta)"  ,Folder),100,-1,1,Folder);
  HBook1F(Hist->fEleMom    ,"ce_mom"   ,Form("%s: Conversion Electron Momentum"    ,Folder),1000,  0,200,Folder);
  HBook1D(Hist->fDioMom    ,"dio_mom"  ,Form("%s: DIO momentum"                    ,Folder),1000, 50,150,Folder);
  HBook1F(Hist->fRv         ,"rv"      ,Form("%s: R(Vertex)"                       ,Folder), 100, 0, 1000,Folder);
  HBook1F(Hist->fZv         ,"zv"      ,Form("%s: Z(Vertex)"                       ,Folder), 300, 0,15000,Folder);
  //  HBook1F(Hist->fNClusters ,"ncl"      ,Form("%s: Number of Reconstructed Clusters",Folder),200,0,200,Folder);
  HBook1F(Hist->fNTracks   ,"ntrk"     ,Form("%s: Number of Reconstructed Tracks"  ,Folder),100,0,100,Folder);
  HBook1F(Hist->fNStrawHits[0],"nsh_0" ,Form("%s: Number of Straw Hits [0]"        ,Folder),250,0,250,Folder);
  HBook1F(Hist->fNStrawHits[1],"nsh_1" ,Form("%s: Number of Straw Hits [1]"        ,Folder),250,0,5000,Folder);
  HBook1F(Hist->fNGoodSH   ,"nsh50"    ,Form("%s: N(SH) +/-50"                     ,Folder),300,0,1500,Folder);
  HBook1F(Hist->fDtClT     ,"dt_clt"   ,Form("%s: DT(cluster-track)"               ,Folder),100,-100,100,Folder);
  HBook1F(Hist->fDtClS     ,"dt_cls"   ,Form("%s: DT(cluster-straw hit)"           ,Folder),200,-200,200,Folder);
  HBook1F(Hist->fSHTime    ,"shtime"   ,Form("%s: Straw Hit Time"                  ,Folder),400,0,2000,Folder);
  HBook1F(Hist->fEMax      ,"emax"     ,Form("%s: Max cluster energy"              ,Folder),150,0,150,Folder);
  HBook1F(Hist->fNHyp      ,"nhyp"     ,Form("%s: N(fit hypotheses)"               ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[0],"bfh0"     ,Form("%s: Best Fit Hyp[0](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fBestHyp[1],"bfh1"     ,Form("%s: Best Fit Hyp[1](e-,e+,mu-,mu+)"  ,Folder),5,0,5,Folder);
  HBook1F(Hist->fNGenp     ,"ngenp"    ,Form("%s: N(Gen Particles)"                ,Folder),500,0,500,Folder);

  // HBook1F(Hist->fNHitCrystalsTot,"nhcr_tot",Form("%s: NHit Crystals Tot",Folder), 100, 0,100,Folder);
  // HBook1F(Hist->fECal,"ecal",Form("%s: E(cal), sum over both disks",Folder), 500, 0,250,Folder);
  // HBook1F(Hist->fECalOverEKin,"ec_over_ek",Form("%s: E(cal)/E(kin)",Folder), 200, 0,2,Folder);
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::BookSimpHistograms(SimpHist_t* Hist, const char* Folder) {
  //  char name [200];
  //  char title[200];

  HBook1F(Hist->fPdgCode   ,"pdg"         ,Form("%s: PDG code"                     ,Folder),200,-100,100,Folder);
  HBook1F(Hist->fNStrawHits,"nsth"        ,Form("%s: n straw hits"                 ,Folder),200,   0,200,Folder);
  HBook1F(Hist->fMomTargetEnd    ,"ptarg" ,Form("%s: CE mom after Stopping Target" ,Folder),400,  90,110,Folder);
  HBook1F(Hist->fMomTrackerFront ,"pfront",Form("%s: CE mom at the Tracker Front"  ,Folder),400,  90,110,Folder);
}


//_____________________________________________________________________________
void TStnPiplusenuAnaModule::BookHistograms() {

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

  book_event_histset[ 0] = 1;  // all events
  book_event_histset[ 1] = 1;  // events with a reconstructed track
  book_event_histset[ 2] = 1;  // events with tan(dip) < 1.1
  book_event_histset[ 3] = 1;  // events with 1.1 < tan(dip) < 1.4
  book_event_histset[ 4] = 1;  // events with tan(dip) > 1.4

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
// book simp histograms
//-----------------------------------------------------------------------------
  int book_simp_histset[kNSimpHistSets];
  for (int i=0; i<kNSimpHistSets; i++) book_simp_histset[i] = 0;

  book_simp_histset[ 0] = 1;		// all events

  for (int i=0; i<kNSimpHistSets; i++) {
    if (book_simp_histset[i] != 0) {
      sprintf(folder_name,"sim_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fSimp[i] = new SimpHist_t;
      BookSimpHistograms(fHist.fSimp[i],Form("Hist/%s",folder_name));
    }
  }
//-----------------------------------------------------------------------------
// book track histograms
//-----------------------------------------------------------------------------
  int book_track_histset[kNTrackHistSets];
  for (int i=0; i<kNTrackHistSets; i++) book_track_histset[i] = 0;

  book_track_histset[  0] = 1;		// all tracks e-
  book_track_histset[  1] = 1;		// all tracks e- passing Set C cuts
  book_track_histset[  2] = 1;		// tracks with tan(dip) < 1.1
  book_track_histset[  3] = 1;		// tracks with 1.1 < tan(dip) < 1.4
  book_track_histset[  4] = 1;		// tracks with tan(dip) > 1.4
  book_track_histset[  5] = 1;		// tracks with tau-weight
  book_track_histset[  6] = 1;		// tracks with t0 < 400 ns

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
// book Genp histograms
//-----------------------------------------------------------------------------
  int book_genp_histset[kNGenpHistSets];
  for (int i=0; i<kNGenpHistSets; i++) book_genp_histset[i] = 0;

  book_genp_histset[0] = 1;		// all particles
  // book_genp_histset[1] = 1;		// all crystals, e > 0
  book_genp_histset[2] = 1;     // tracks with tan(dip) < 1.1
  book_genp_histset[3] = 1;     // tracks with 1.1 < tan(dip) < 1.4
  book_genp_histset[4] = 1;     // tracks with tan(dip) > 1.4

  for (int i=0; i<kNGenpHistSets; i++) {
    if (book_genp_histset[i] != 0) {
      sprintf(folder_name,"gen_%i",i);
      fol = (TFolder*) hist_folder->FindObject(folder_name);
      if (! fol) fol = hist_folder->AddFolder(folder_name,folder_name);
      fHist.fGenp[i] = new GenpHist_t;
      BookGenpHistograms(fHist.fGenp[i],Form("Hist/%s",folder_name));
    }
  }
}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::FillEventHistograms(EventHist_t* Hist) {
  double            cos_th(-2), dio_wt(-1.), xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6), p(-1.);
  //  double            e, m, r;
  TLorentzVector    mom;

  if (fParticle) {
    fParticle->Momentum(mom);
    p      = mom.P();
    cos_th = mom.Pz()/p;
    xv     = fParticle->Vx()+3904.;
    yv     = fParticle->Vy();
    rv     = sqrt(xv*xv+yv*yv);
    zv     = fParticle->Vz();
    dio_wt = TStntuple::DioWeightAl(p);
  }

  Hist->fEleMom->Fill(p);
  Hist->fDioMom->Fill(p,dio_wt);
  Hist->fEleCosTh->Fill(cos_th);
  Hist->fRv->Fill(rv);
  Hist->fZv->Fill(zv);

  // Hist->fNClusters->Fill(fNClusters);
  Hist->fNTracks->Fill  (fNTracks[0]);
  // Hist->fNStrawHits[0]->Fill(fNStrawHits);
  // Hist->fNStrawHits[1]->Fill(fNStrawHits);

  double emax   = -1;
  double dt     = 9999.;

  Hist->fDtClT->Fill(dt);
  Hist->fEMax->Fill(emax);

  Hist->fNHyp->Fill(fNHyp);
  Hist->fBestHyp[0]->Fill(fBestHyp[0]);
  Hist->fBestHyp[1]->Fill(fBestHyp[1]);
  Hist->fNGenp->Fill(fNGenp);
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::FillGenpHistograms(GenpHist_t* Hist, TGenParticle* Genp) {
  int    gen_id;
  float  p, cos_th, z0, t0, r0, x0, y0;

  TLorentzVector mom, v;

  Genp->Momentum(mom);
  //  Genp->ProductionVertex(v);

  p      = mom.P();
  cos_th = mom.CosTheta();

  x0     = Genp->Vx()+3904.;
  y0     = Genp->Vy();

  z0     = Genp->Vz();
  t0     = Genp->T();
  r0     = sqrt(x0*x0+y0*y0);
  gen_id = Genp->GetStatusCode();

  Hist->fPdgCode[0]->Fill(Genp->GetPdgCode());
  Hist->fPdgCode[1]->Fill(Genp->GetPdgCode());
  Hist->fGenID->Fill(gen_id);
  Hist->fX0->Fill(Genp->Vx());
  Hist->fY0->Fill(y0);
  Hist->fZ0->Fill(z0);
  Hist->fT0->Fill(t0);
  Hist->fR0->Fill(r0);
  Hist->fP->Fill(p);
  Hist->fCosTh->Fill(cos_th);
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::FillSimpHistograms(SimpHist_t* Hist, TSimParticle* Simp) {

  Hist->fPdgCode->Fill(Simp->fPdgCode);
  Hist->fMomTargetEnd->Fill(Simp->fMomTargetEnd);
  Hist->fMomTrackerFront->Fill(Simp->fMomTrackerFront);
  Hist->fNStrawHits->Fill(Simp->fNStrawHits);
}

//-----------------------------------------------------------------------------
// for DIO : ultimately, one would need to renormalize the distribution
//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::FillTrackHistograms(TrackHist_t* Hist, TStnTrack* Track, double Weight) {

  TLorentzVector  mom;
  double          r;
  int             itrk;
  TrackPar_t*     tp;
					// pointer to local track parameters
  itrk = Track->Number();
  tp   = fTrackPar+itrk;

  Hist->fP[0]->Fill (Track->fP);
  Hist->fP[1]->Fill (Track->fP);
  Hist->fP[2]->Fill (Track->fP);
  Hist->fP0->  Fill (Track->fP0);
  Hist->fP2->  Fill (Track->fP2);

  Hist->fPDio->Fill(Track->fP,Weight);

  Hist->fFitMomErr->Fill(Track->fFitMomErr,Weight);

  Hist->fPt    ->Fill(Track->fPt    ,Weight);
  Hist->fPFront->Fill(Track->fPFront,Weight);
  Hist->fPStOut->Fill(Track->fPStOut,Weight);
					// dp: Tracker-only resolution
  Hist->fDpFront ->Fill(tp->fDpF,Weight);
  Hist->fDpFront0->Fill(tp->fDp0,Weight);
  Hist->fDpFront2->Fill(tp->fDp2,Weight);
  Hist->fDpFSt   ->Fill(tp->fDpFSt,Weight);
  Hist->fDpFVsZ1 ->Fill(Track->fZ1,tp->fDpF,Weight);

  Hist->fCosTh->Fill(Track->Momentum()->CosTheta(),Weight);
  Hist->fChi2->Fill (Track->fChi2,Weight);
  Hist->fNDof->Fill(Track->NActive()-5.,Weight);
  Hist->fChi2Dof->Fill(Track->fChi2/(Track->NActive()-5.),Weight);
  Hist->fNActive->Fill(Track->NActive(),Weight);
  Hist->fT0->Fill(Track->fT0,Weight);
  Hist->fT0Err->Fill(Track->fT0Err,Weight);
  Hist->fQ->Fill(Track->Charge(),Weight);
  Hist->fFitCons[0]->Fill(Track->fFitCons,Weight);
  Hist->fFitCons[1]->Fill(Track->fFitCons,Weight);

  Hist->fD0->Fill(Track->fD0,Weight);
  Hist->fZ0->Fill(Track->fZ0,Weight);
  Hist->fTanDip->Fill(Track->fTanDip,Weight);
  Hist->fAlgMask->Fill(Track->AlgMask(),Weight);

//-----------------------------------------------------------------------------
// track-cluster matching part: 
// - for residuals, determine intersection with the most energetic cluster
// - for track -only parameters use intersection with lowest trajectory length
//-----------------------------------------------------------------------------
  TStnTrack::InterData_t*    vt = Track->fVMinS;  // track-only

  if (vt) {
    Hist->fVaneID->Fill(vt->fID  ,Weight);
    Hist->fXTrk->Fill  (vt->fXTrk,Weight);
    Hist->fYTrk->Fill  (vt->fYTrk,Weight);

    r = sqrt(vt->fXTrk*vt->fXTrk+vt->fYTrk*vt->fYTrk);
    Hist->fRTrk->Fill  (r,Weight);

    Hist->fZTrk->Fill  (vt->fZTrk,Weight);
  }
  else {
//-----------------------------------------------------------------------------
// fill histograms with numbers easy to recognize as dummy
//-----------------------------------------------------------------------------
    Hist->fVaneID->Fill( -1.,Weight);
    Hist->fXTrk->Fill  (999.,Weight);
    Hist->fYTrk->Fill  (999.,Weight);
    Hist->fRTrk->Fill  (999.,Weight);
    Hist->fZTrk->Fill  ( -1.,Weight);
  }
//-----------------------------------------------------------------------------
// there is an inconsistency in the SIMP block filling - in Mu2e offline 
// the particle momentumis is kept in MeV/c, while the PDG mass  -in GeV/c^2..
// thus the energy is screwed up... kludge around
// assign muon mass
//-----------------------------------------------------------------------------
  double ekin(-1.);
  if (fSimp) {
    double p, m;
    //    p    = fSimp->fStartMom.P();
    p = Track->fP;
    m    = 105.658; // in MeV
    ekin = sqrt(p*p+m*m)-m;
  }

  Hist->fECl->Fill(tp->fEcl,Weight);
  Hist->fEClEKin->Fill(tp->fEcl/ekin,Weight);
  Hist->fEp->Fill(tp->fEp,Weight);
  Hist->fEpVsPath->Fill(tp->fPath,tp->fEp,Weight);

  Hist->fDx->Fill(tp->fDx,Weight);
  Hist->fDy->Fill(tp->fDy,Weight);
  Hist->fDz->Fill(tp->fDz,Weight);

  Hist->fDt->Fill(tp->fDt,Weight);
  Hist->fChi2Match->Fill(tp->fChi2Match,Weight);

  Hist->fDu->Fill    (tp->fDu,Weight);
  Hist->fDv->Fill    (tp->fDv,Weight);
  Hist->fDvVsDu->Fill(tp->fDu,tp->fDv,Weight);

  Hist->fPath->Fill(tp->fPath,Weight);
  Hist->fDuVsPath->Fill(tp->fPath,tp->fDu,Weight);

  //  double cu[4] = { -60.3049, -0.749111,    0.00522242,  -7.52018e-06};
  double cu[4] = { -59.5174, -0.541226, 0.00414309, -5.84989e-06 };
  double cv[4] = {  6.44161, 0.0722353, -0.000653084, 1.14054e-06};

  double x = tp->fPath;
  double corr_u = cu[0]+cu[1]*x+cu[2]*x*x+cu[3]*x*x*x;
  double corr_v = cv[0]+cv[1]*x+cv[2]*x*x+cv[3]*x*x*x;

  //  double duc = tp->fDu-0.34*(tp->fPath-350.);
  double duc = tp->fDu-corr_u;
  double dvc = tp->fDv-corr_v;
  
  Hist->fDucVsPath->Fill(tp->fPath,duc,Weight);
  Hist->fDvcVsPath->Fill(tp->fPath,dvc,Weight);

  Hist->fDvVsPath->Fill(tp->fPath,tp->fDv,Weight);
  Hist->fDtVsPath->Fill(tp->fPath,tp->fDt,Weight);

  Hist->fDuVsTDip->Fill(Track->fTanDip,tp->fDu,Weight);
  Hist->fDvVsTDip->Fill(Track->fTanDip,tp->fDv,Weight);

  Hist->fZ1->Fill(Track->fZ1);

  int ncl = Track->NClusters();
  Hist->fNClusters->Fill(ncl);

  Hist->fRSlope->Fill(Track->RSlope(),Weight);
  Hist->fXSlope->Fill(Track->XSlope(),Weight);

  double llhr_dedx, llhr_xs, llhr_cal, llhr_trk, llhr;

  Hist->fEleLogLHCal->Fill(Track->EleLogLHCal(),Weight);
  Hist->fMuoLogLHCal->Fill(Track->MuoLogLHCal(),Weight);

  llhr_cal = Track->LogLHRCal();
  Hist->fLogLHRCal->Fill(llhr_cal,Weight);

  llhr_dedx = Track->LogLHRDeDx();
  llhr_xs   = Track->LogLHRXs();
  llhr_trk  = Track->LogLHRTrk();
  llhr      = llhr_cal+llhr_trk;

  Hist->fEpVsDt->Fill(tp->fDt,tp->fEp,Weight);
  Hist->fLogLHRDeDx->Fill(llhr_dedx,Weight);
  Hist->fLogLHRXs->Fill(llhr_xs,Weight);
  Hist->fLogLHRTrk->Fill(llhr_trk,Weight);
  Hist->fLogLHR->Fill(llhr,Weight);

  Hist->fPdgCode->Fill(Track->fPdgCode,Weight);
  Hist->fFrGH->Fill(Track->fNGoodMcHits/(Track->NActive()+1.e-5),Weight);

  Hist->fNEPlVsNHPl->Fill(tp->fNEPl,tp->fNHPl,Weight);
  Hist->fNDPlVsNHPl->Fill(tp->fNDPl,tp->fNHPl,Weight);
  Hist->fChi2dVsNDPl->Fill(tp->fNDPl,Track->Chi2Dof(),Weight);
  Hist->fDpFVsNDPl  ->Fill(tp->fNDPl,tp->fDpF,Weight);
}


//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TStnPiplusenuAnaModule::BeginJob() {
//-----------------------------------------------------------------------------

// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("TrackBlockDar"       ,"TStnTrackBlock"      ,&fTrackBlock      );
  RegisterDataBlock("GenpBlock"           ,"TGenpBlock"          ,&fGenpBlock       );
  RegisterDataBlock("SimpBlock"           ,"TSimpBlock"          ,&fSimpBlock       );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms();
//-----------------------------------------------------------------------------
// initialize likelihood histograms
//-----------------------------------------------------------------------------
  fTrackID->SetMinT0(fMinT0);
//-----------------------------------------------------------------------------
// PID initialization: read the likelihood templates
//-----------------------------------------------------------------------------
  fLogLH->Init("v5_7_0");

  return 0;
}


//_____________________________________________________________________________
int TStnPiplusenuAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TStnPiplusenuAnaModule::FillHistograms() {

  //  double       cos_th (-2.),  cl_e(-1.);
  // int          disk_id(-1); // , alg_mask, nsh, nactive;
  //  float        pfront, ce_pitch, reco_pitch, fcons, t0, sigt, sigp, p; 
  //  TStnCluster  *cl0;

  //  cos_th = fEle->momentum().pz()/fEle->momentum().vect().mag();

  // if (fNClusters > 0) {
  //   cl0     = fClusterBlock->Cluster(0);
  //   cl_e    = cl0->Energy();
  //   disk_id = cl0->DiskID();
  // }
  //-----------------------------------------------------------------------------
  // event histograms
  //-----------------------------------------------------------------------------
  FillEventHistograms(fHist.fEvent[0]);

  if (fNTracks[0]> 0) FillEventHistograms(fHist.fEvent[1]);
//-----------------------------------------------------------------------------
// Simp histograms
//-----------------------------------------------------------------------------
  if (fSimp) {
    FillSimpHistograms(fHist.fSimp[0],fSimp);
  }
//-----------------------------------------------------------------------------
// track histograms, fill them only for the downstream e- hypothesis
//-----------------------------------------------------------------------------
  TStnTrack*   trk;
  //  TrackPar_t*  tp;
 
  for (int i=0; i<fNTracks[0]; ++i ) {
    trk = fTrackBlock->Track(i);
    //    tp  = fTrackPar+i;

    FillTrackHistograms(fHist.fTrack[0],trk);

    if (trk->fIDWord == 0) {
    // track passes selection "C"

      FillTrackHistograms(fHist.fTrack[1],trk);
    }

    if (trk->fTanDip < 1.1) {
    }

    if (trk->T0() < 400){
      FillTrackHistograms(fHist.fTrack[6],trk);
    }
  }
//-----------------------------------------------------------------------------
// fill GENP histograms
// GEN_0: all particles
//-----------------------------------------------------------------------------
  TGenParticle* genp;
  for (int i=0; i<fNGenp; i++) {
    genp = fGenpBlock->Particle(i);
    FillGenpHistograms(fHist.fGenp[0],genp);

    // Fill Genp, Track, and Event histograms based on selections
    if (fNTracks[0]> 0) {
      for (int k=0; k<fNTracks[0]; ++k) {
        trk = fTrackBlock->Track(k);
        double weight = fGenpBlock->Weight();
        FillTrackHistograms(fHist.fTrack[5],trk,weight);

        if (trk->fTanDip < 1.1) {
          FillGenpHistograms(fHist.fGenp[2],genp);
          FillTrackHistograms(fHist.fTrack[2],trk);
          FillEventHistograms(fHist.fEvent[2]);
        }
        else if (trk->fTanDip > 1.1 && trk->fTanDip < 1.4) {
          FillGenpHistograms(fHist.fGenp[3],genp);
          FillTrackHistograms(fHist.fTrack[3],trk);
          FillEventHistograms(fHist.fEvent[3]);
        }
        else if (trk->fTanDip > 1.4) {
          FillGenpHistograms(fHist.fGenp[4],genp);
          FillTrackHistograms(fHist.fTrack[4],trk);
          FillEventHistograms(fHist.fEvent[4]);
        }
      }
    }
  }
  //  first_entry = 0;
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TStnPiplusenuAnaModule::Event(int ientry) {

  double                xs, p;
  TEmuLogLH::PidData_t  dat;
  TStnTrack*            track;
  int                   id_word;
  TLorentzVector        mom;

  //  TDiskCalorimeter::GeomData_t disk_geom;

  fTrackBlock  ->GetEntry(ientry);
  fGenpBlock->GetEntry(ientry);
  fSimpBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNGenp    = fGenpBlock->NParticles();

  TGenParticle* genp;
  int           pdg_code, generator_code;

  fParticle = NULL;
  for (int i=fNGenp-1; i>=0; i--) {
    genp           = fGenpBlock->Particle(i);
    pdg_code       = genp->GetPdgCode();
    generator_code = genp->GetStatusCode();
    if ((abs(pdg_code) == fPdgCode) && (generator_code == fGeneratorCode)) {
      fParticle = genp;
      break;
    }
  }
					// may want to revisit the definition of fSimp
  fSimp     = fSimpBlock->Particle(0);

  if (fParticle) {
    fParticle->Momentum(mom);
    p         = mom.P();
  }
  else p = 0.;

  fEleE     = sqrt(p*p+0.511*0.511);


  fNTracks[0] = fTrackBlock->NTracks();

  fNHyp       = -1;
  fBestHyp[0] = -1;
  fBestHyp[1] = -1;

  fNGoodTracks    = 0;
  fNMatchedTracks = 0;

  fNTracks[0] = fTrackBlock->NTracks();
  if (fNTracks[0] == 0) fTrack = 0;
  else                  fTrack = fTrackBlock->Track(0);

  int ntrk = fNTracks[0];

  TrackPar_t*   tp;

  for (int itrk=0; itrk<ntrk; itrk++) {
					// assume less 20 tracks
    tp             = fTrackPar+itrk;

    track          = fTrackBlock->Track(itrk);
    id_word        = fTrackID->IDWord(track);
    track->fIDWord = id_word;
    if (id_word == 0) {
      fNGoodTracks += 1;
      if ((track->fVMaxEp != NULL) && (fabs(track->fVMaxEp->fDt) < 2.5)) {
	fNMatchedTracks += 1;
      }
    }
//-----------------------------------------------------------------------------
// process hit masks
//-----------------------------------------------------------------------------
    int i1, i2, n1(0) ,n2(0), ndiff(0);
    int nbits = track->fHitMask.GetNBits();
    for (int i=0; i<nbits; i++) {
      i1 = track->HitMask()->GetBit(i);
      i2 = track->ExpectedHitMask()->GetBit(i);
      n1 += i1;
      n2 += i2;
      if (i1 != i2) ndiff += 1;
    }
//-----------------------------------------------------------------------------
// define additional parameters
//-----------------------------------------------------------------------------
    tp->fNHPl = n1;
    tp->fNEPl = n2;
    tp->fNDPl = ndiff;

    tp->fDpF   = track->fP     -track->fPFront;
    tp->fDp0   = track->fP0    -track->fPFront;
    tp->fDp2   = track->fP2    -track->fPFront;
    tp->fDpFSt = track->fPFront-track->fPStOut;

    if (fFillDioHist == 0) tp->fDioWt = 1.;
    else                   tp->fDioWt = TStntuple::DioWeightAl(fEleE);
//-----------------------------------------------------------------------------
// track residuals
//-----------------------------------------------------------------------------
    TStnTrack::InterData_t*  vr = track->fVMaxEp; 
    double    nx, ny;

    tp->fEcl       = -1.e6;
    tp->fEp        = -1.e6;

    tp->fDu        = -1.e6;
    tp->fDv        = -1.e6;
    tp->fDx        = -1.e6;
    tp->fDy        = -1.e6;
    tp->fDz        = -1.e6;
    tp->fDt        = -1.e6;

    tp->fChi2Match = -1.e6;
    tp->fPath      = -1.e6;

    if (vr) {
      tp->fEcl = vr->fEnergy;
      tp->fEp  = tp->fEcl/track->fP;

      tp->fDx  = vr->fDx;
      tp->fDy  = vr->fDy;
      tp->fDz  = vr->fDz;
//-----------------------------------------------------------------------------
// v4_2_4: correct by additional 0.22 ns - track propagation by 6 cm
//-----------------------------------------------------------------------------
      tp->fDt  = vr->fDt - 0.22; // - 1.;

      nx  = vr->fNxTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);
      ny  = vr->fNyTrk/sqrt(vr->fNxTrk*vr->fNxTrk+vr->fNyTrk*vr->fNyTrk);

      tp->fDu        = vr->fDx*nx+vr->fDy*ny;
      tp->fDv        = vr->fDx*ny-vr->fDy*nx;
      tp->fChi2Match = vr->fChi2Match;
      tp->fPath      = vr->fPath;
    }

    //    if ((tp->fEp > 0) && (track->fEp > 0) && (fabs(tp->fEp-track->fEp) > 1.e-6)) {
    //      GetHeaderBlock()->Print(Form(" TStnPiplusenuAnaModule ERROR: tp->fEp = %10.5f  track->fEp = %10.5f",tp->fEp,track->fEp));
    //    }
//-----------------------------------------------------------------------------
// PID likelihoods
//-----------------------------------------------------------------------------
    dat.fDt   = tp->fDt;
    dat.fEp   = tp->fEp;
    dat.fPath = tp->fPath;
      
    xs = track->XSlope();

    track->fEleLogLHCal = fLogLH->LogLHCal(&dat,11);
    track->fMuoLogLHCal = fLogLH->LogLHCal(&dat,13);

    double llhr_cal = track->fEleLogLHCal-track->fMuoLogLHCal;

    if (GetDebugBit(7)) {
      if ((id_word == 0) && (llhr_cal > 20)) {
	GetHeaderBlock()->Print(Form("bit:007: dt = %10.3f ep = %10.3f",track->Dt(),tp->fEp));
      }
    }

    if (GetDebugBit(8)) {
      if ((id_word == 0) && (llhr_cal < -20)) {
	GetHeaderBlock()->Print(Form("bit:008: p = %10.3f dt = %10.3f ep = %10.3f",
				     track->P(),track->Dt(),tp->fEp));
      }
    }

    track->fLogLHRXs    = fLogLH->LogLHRXs(xs);
  }

  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TStnPiplusenuAnaModule::Debug() {

  TStnTrack* trk;
  //  TrackPar_t* tp;
  int ntrk = fTrackBlock->NTracks();

  for (int itrk=0; itrk<ntrk; itrk++) {
    trk = fTrackBlock->Track(itrk);
    //    tp  = &fTrackPar[itrk];
//-----------------------------------------------------------------------------
// bit 3: Set C tracks with large DX : 70mm < |DX| < 90mm
//-----------------------------------------------------------------------------
    if (GetDebugBit(3) == 1) {
      if (trk->fIDWord == 0) {
	TStnTrack::InterData_t*    vr = trk->fVMaxEp; // residuals
	if ((vr && (fabs(vr->fDx) > 70) && (fabs(vr->fDx) < 90))) {
	  GetHeaderBlock()->Print(Form("large DX: %f",vr->fDx));
	}
      }
    }

    if (GetDebugBit(4) == 1) {
      if (trk->TanDip() > 1.4) {
	GetHeaderBlock()->Print(Form("large tan(dip): %f",trk->TanDip()));
      }
    }
    if (GetDebugBit(5) == 1) {
      if (trk->P() > 60) {
	GetHeaderBlock()->Print(Form(" trk p > 60: %f",trk->P()));
      }
    }
  }
}

//_____________________________________________________________________________
int TStnPiplusenuAnaModule::EndJob() {
  printf("----- end job: ---- %s\n",GetName());
  return 0;
}

//_____________________________________________________________________________
void TStnPiplusenuAnaModule::Test001() {
}
