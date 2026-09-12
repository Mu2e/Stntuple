//////////////////////////////////////////////////////////////////////////////
// use of tmp:
//
// use of debug bits: bits 0-2 are reserved
//  0  : all events
//  1  : passed events
//  2  : rejected events
//  3  : N(CalHelixFinder helices hel>0) > 0
//  4  : N(CalHelixFinder helices hel<0) > 0
// 
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <format>

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
#include "ana/TDetTimeAnaModule.hh"


ClassImp(stntuple::TDetTimeAnaModule)

namespace stntuple {
//-----------------------------------------------------------------------------
TDetTimeAnaModule::TDetTimeAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  TH1::AddDirectory(0);
//-----------------------------------------------------------------------------
// detector numerology
//-----------------------------------------------------------------------------
  fTpm            = TrkPanelMap::Instance();
  fCaloChannelMap = TCaloChannelMap::Instance();
  
  fFrRef          = nullptr;
  fHist           = new Hist_t;
}

//-----------------------------------------------------------------------------
TDetTimeAnaModule::~TDetTimeAnaModule() {
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::BookCalcHistograms(CalcHist_t* Hist, CaloIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  name  = "edep";
  title = std::format("{} : edep",prefix);
  fBookHist->HBook1F(Hist->h_edep,name.data(),title.data(),100,0,1000,Folder);   // in MeV

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBookHist->HBook1F(Hist->h_dt_tc,name.data(),title.data(),100,-100,100,Folder);   // in ns...

  name  = "dt_crvc";
  title = std::format("{} : dt CRVC",prefix);
  fBookHist->HBook1F(Hist->h_dt_crvc,name.data(),title.data(),100,-100,100,Folder);   // in ns...

  name  = "dt_crvc_vs_dt_tc";
  title = std::format("{} : dt CRVC vs dt TC",prefix);
  fBookHist->HBook2F(Hist->h_dt_crvc_vs_dt_tc,name.data(),title.data(),100,-100,100,100,-100,100,Folder);   // in ns...

  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::BookCalhHistograms(CalhHist_t* Hist, CaloIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} disk:{:02d} crate:{}",
                                   fRunNumber,Index->sel, Index->disk, Index->crate);
  std::string name, title;

  // name  = "ph";
  // title = std::format("{} : ph",prefix);
  // fBookHist->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...



  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::BookTrkHistograms(TrkHist_t* Hist, TrkIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "nhits";
  title = std::format("{} : nhits",prefix);
  fBookHist->HBook1F(Hist->h_nhits,name.data(),title.data(),100,0,100,Folder);   // in MeV

  name  = "chi2d";
  title = std::format("{} : chi2/ndod",prefix);
  fBookHist->HBook1F(Hist->h_chi2d,name.data(),title.data(),100,0,20,Folder);   // in ns...

  name  = "t0";
  title = std::format("{} : t0",prefix);
  fBookHist->HBook1F(Hist->h_t0,name.data(),title.data(),500,0,2.5e6,Folder);   // in ns...

  name  = "dt_tc";
  title = std::format("{} : dt TC",prefix);
  fBookHist->HBook1F(Hist->h_dt_tc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dt_calc";
  title = std::format("{} : dt CALC",prefix);
  fBookHist->HBook1F(Hist->h_dt_calc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dt_crvc";
  title = std::format("{} : dt CRVC",prefix);
  fBookHist->HBook1F(Hist->h_dt_crvc,name.data(),title.data(),400,-100,100,Folder);   // in ns...

  name  = "dx_calc";
  title = std::format("{} : dx_calc = x(trk)-x(calc)",prefix);
  fBookHist->HBook1F(Hist->h_dx_calc,name.data(),title.data(),200,-1000,1000,Folder);   // in ns...

  name  = "dy_calc";
  title = std::format("{} : dy_calc = y(trk)-y(calc)",prefix);
  fBookHist->HBook1F(Hist->h_dy_calc,name.data(),title.data(),100,-2500,2500,Folder);   // in ns...

  name  = "dx_calc_vs_dxdz";
  title = std::format("{} : dx_calc vs dxdz",prefix);
  fBookHist->HBook2F(Hist->h_dx_calc_vs_dxdz,name.data(),title.data(),200,-1,1,200,-1000,1000,Folder);

  name  = "dy_calc_vs_dydz";
  title = std::format("{} : dy_calc vs dydz",prefix);
  fBookHist->HBook2F(Hist->h_dy_calc_vs_dydz,name.data(),title.data(),100,-2.5,2.5,200,-1000,1000,Folder);

  name  = "dxdz";
  title = std::format("{} : dxdz",prefix);
  fBookHist->HBook1F(Hist->h_dxdz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "dydz";
  title = std::format("{} : dydz",prefix);
  fBookHist->HBook1F(Hist->h_dydz,name.data(),title.data(),200,-10,10,Folder);   // in ns...

  name  = "xcrv";
  title = std::format("{} : xcrv",prefix);
  fBookHist->HBook1F(Hist->h_xcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "zcrv";
  title = std::format("{} : zcrv",prefix);
  fBookHist->HBook1F(Hist->h_zcrv,name.data(),title.data(),200,-10000,10000,Folder);   // in ns...

  name  = "dx_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBookHist->HBook1F(Hist->h_dx_crvc,name.data(),title.data(),200,-10000,10000,Folder);   //

  name  = "dz_crvc";
  title = std::format("{} : dx_crvc",prefix);
  fBookHist->HBook1F(Hist->h_dz_crvc,name.data(),title.data(),200,-10000,10000,Folder);   //

  name  = "dx_crvc_vs_dxdy";
  title = std::format("{} : dx_crvc vs dxdy",prefix);
  fBookHist->HBook2F(Hist->h_dx_crvc_vs_dxdy,name.data(),title.data(),100,-2,2,200,-1000,1000,Folder);

  name  = "dz_crvc_vs_dzdy";
  title = std::format("{} : dz_crvc vs dzdy",prefix);
  fBookHist->HBook2F(Hist->h_dz_crvc_vs_dzdy,name.data(),title.data(),100,-2,2,200,-1000,1000,Folder);

  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CaloIndex_t index;

  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "dt_vs_sipmid";
  title = std::format("{} : dt vs crystal ID",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_sipmid,name.data(),title.data(),3000,0,3000,1000,-1000,1000,Folder);

  name  = "sipmid";
  title = std::format("{} : Sipm ID",prefix);
  fBookHist->HBook1F(Hist->h_sipmid,name.data(),title.data(),3000,0,3000,Folder);

  name  = "n2_vs_n1";
  title = std::format("{} : N1:N1 calh E>10",prefix);
  fBookHist->HBook2F(Hist->h_n2_vs_n1,name.data(),title.data(),20,0,20,20,0,20,Folder);

  name  = "ntrk_0";
  title = std::format("{} : ntrk[0]",prefix);
  fBookHist->HBook1F(Hist->h_ntrk[0],name.data(),title.data(),10,0,10,Folder);

  name  = "ntrk_1";
  title = std::format("{} : ntrk[1]",prefix);
  fBookHist->HBook1F(Hist->h_ntrk[1],name.data(),title.data(),10,0,10,Folder);
//-----------------------------------------------------------------------------
// calorimeter hits 
//-----------------------------------------------------------------------------
  int book_calh_histset[kNCalhHistSets];

  for (int i=0; i<kNCalhHistSets; i++) { book_calh_histset[i] = 0; }

  book_calh_histset[0] = 1;             // all

  for (int i=0; i<kNCalhHistSets; i++) {
    if (book_calh_histset[i] == 0) continue;
    std::string folder_name = std::format("calh_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calh[i] = new CalhHist_t();
    index.sel = i;
    BookCalhHistograms(Hist->calh[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CALC: calorimeter clusters 
//-----------------------------------------------------------------------------
  int book_calc_histset[kNCalcHistSets];
  for (int i=0; i<kNCalcHistSets; i++) { book_calc_histset[i] = 0; }

  book_calc_histset[0] = 1;             // all
  book_calc_histset[1] = 1;             // all DISK 0
  book_calc_histset[2] = 1;             // all DISK 1
  book_calc_histset[3] = 1;             // |dt_tc| < 30, all 
  book_calc_histset[4] = 1;             // |dt_tc| < 30, DISK 0
  book_calc_histset[5] = 1;             // |dt_tc| < 30, DISK 1
  book_calc_histset[6] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, all
  book_calc_histset[7] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, DISK 0
  book_calc_histset[8] = 1;             // |dt_tc| < 30, |dt_crvc| < 30, DISK 1

  int sel_disk[100];
  sel_disk[0] = -1;
  sel_disk[1] =  0;
  sel_disk[2] =  1;
  sel_disk[3] = -1;
  sel_disk[4] =  0;
  sel_disk[5] =  1;
  sel_disk[6] = -1;
  sel_disk[7] =  0;
  sel_disk[8] =  1;

  for (int i=0; i<kNCalcHistSets; i++) {
    if (book_calc_histset[i] == 0) continue;
    std::string folder_name = std::format("calc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->calc[i] = new CalcHist_t();
    index.sel  = i;
    index.disk = sel_disk[i];
    BookCalcHistograms(Hist->calc[i],&index,fol);
  }
//-----------------------------------------------------------------------------
// Trk: tracks 
//-----------------------------------------------------------------------------
  TrkIndex_t trk_index;

  // std::string prefix = std::format("run:{:06d}",fRunNumber);
  int book_trk_histset[kNTrkHistSets];

  for (int i=0; i<kNTrkHistSets; i++) { book_trk_histset[i] = 0; }

  book_trk_histset[0] = 1;             // all
  
  book_trk_histset[1] = 1;             // in-time tracks
  book_trk_histset[2] = 1;             // in-time tracks DISK0
  book_trk_histset[3] = 1;             // in-time tracks DISK1

  book_trk_histset[4] = 1;             // in-time CAL
  book_trk_histset[5] = 1;             // in-time CAL DISK0
  book_trk_histset[6] = 1;             // in-time CAL DISK1

  book_trk_histset[7] = 1;             // in-time trk_nhits>10 CAL
  book_trk_histset[8] = 1;             // in-time trk_nhits>10 CAL DISK0
  book_trk_histset[9] = 1;             // in-time trk_nhits>10 CAL DISK1

  book_trk_histset[10] = 1;             // in-time CRVC trk_nhits>10

  for (int i=0; i<kNTrkHistSets; i++) {
    if (book_trk_histset[i] == 0) continue;
    std::string folder_name = std::format("trk_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->trk[i]   = new TrkHist_t();
    trk_index.sel  = i;
    BookTrkHistograms(Hist->trk[i],&trk_index,fol);
  }

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TDetTimeAnaModule::FillCalcHistograms(CalcHist_t* Hist, TStnCluster* Calc, calc_param_t* Cp) {
  // filling histograms: plot time differences between
  
  Hist->h_edep->Fill(Calc->Energy());
  Hist->h_dt_tc->Fill(Cp->dtmin_tc);
  Hist->h_dt_crvc->Fill(Cp->dtmin_crvc);
  Hist->h_dt_crvc_vs_dt_tc->Fill(Cp->dtmin_tc,Cp->dtmin_crvc);

  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TDetTimeAnaModule::FillDiskHistograms(DiskHist_t* Hist, TCaloRecoDigi* Calrd) {
  // filling histograms: plot time differences between
  // Hist->h_ph->Fill(Crvp->ph);
  // Hist->h_npes->Fill(Crvp->npes);
  // Hist->h_time->Fill(Crvp->time);
  // Hist->h_feb->Fill(Crvp->feb);
  // Hist->h_ch->Fill(Crvp->ch);
  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::FillTrkHistograms(TrkHist_t* Hist, TStrTrack* Trk, trk_param_t* Tp) {
  // filling histograms: plot time differences between
  
  Hist->h_nhits->Fill(Trk->fNHits);
  Hist->h_chi2d->Fill(Trk->fChi2/Trk->fNDof);
  Hist->h_t0->Fill(Trk->fT0);
  Hist->h_dt_tc->Fill(Tp->dtmin_tc);

  Hist->h_dt_calc->Fill(Tp->dtmin_calc);
  Hist->h_dt_crvc->Fill(Tp->dtmin_crvc);

  Hist->h_dx_calc->Fill(Tp->dx_calc);
  Hist->h_dy_calc->Fill(Tp->dy_calc);

  float dxdz = Trk->fNx/Trk->fNz;
  float dydz = Trk->fNy/Trk->fNz;
  
  Hist->h_dx_calc_vs_dxdz->Fill(dxdz,Tp->dx_calc);
  Hist->h_dy_calc_vs_dydz->Fill(dydz,Tp->dy_calc);

  float dxdy = dxdz/dydz;
  float dzdy = 1./dydz;

  Hist->h_dxdz->Fill(dxdz);
  Hist->h_dydz->Fill(dydz);
  
  Hist->h_xcrv->Fill(Tp->xcrv);
  Hist->h_zcrv->Fill(Tp->zcrv);

  Hist->h_dx_crvc->Fill(Tp->dx_crvc);
  Hist->h_dz_crvc->Fill(Tp->dz_crvc);

  if (Tp->crvc) {
    int s_type = Tp->crvc->SectorType();
    if ((s_type == 1) or (s_type == 3)) {
      // EX or T2: long side of bars along X axis, measuring Z
      Hist->h_dz_crvc_vs_dzdy->Fill(dzdy,Tp->dz_crvc);
    }
    else {
      // T1 or M1-M8 : long side of bars along Z axis, the bars tell X
      Hist->h_dx_crvc_vs_dxdy->Fill(dxdy,Tp->dx_crvc);
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
int TDetTimeAnaModule::FillHistograms() {
  // filling histograms: plot time differences between

  //  Index_t index;
  
//-----------------------------------------------------------------------------
// double-nested loops start here
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNCalrd; i1++) {
    TCaloRecoDigi*  calrd = fCaloRecoDigiBlock->CaloRecoDigi(i1);

    // find closest time cluster
    float dtmin = 1.e6;
    for (int i2=0; i2<fNTc; i2++) {
      TStnTimeCluster*  tc = fTcBlock->TimeCluster(i2);
      if (tc->NStrawHits() < 8) continue;

      float dt       = calrd->Time()-tc->T0();
      if (fabs(dt) < dtmin) {
        dtmin = dt;
      }
    }
    
    fHist->h_dt_vs_sipmid->Fill(calrd->SipmID(),dtmin);
    fHist->h_sipmid->Fill(calrd->SipmID());
  }
  
  fHist->h_n2_vs_n1->Fill(fNCalh10[0],fNCalh10[1]);

//-----------------------------------------------------------------------------
// fill cluster histograms
//-----------------------------------------------------------------------------
  int n_good_tc = 0;
  
  for (int i=0; i<fNCaloClusters; i++) {
    TStnCluster* calc        = fCaloClusterBlock->Cluster(i);
    calc_param_t*   calc_par = &fListOfCalcParam[i];
    FillCalcHistograms(fHist->calc[0],calc,calc_par);
    if (calc->DiskID() == 0) FillCalcHistograms(fHist->calc[1],calc,calc_par);
    else                     FillCalcHistograms(fHist->calc[2],calc,calc_par);

    if (fabs(calc_par->dtmin_tc) < 30) {
      FillCalcHistograms(fHist->calc[3],calc,calc_par);
      if (calc->DiskID() == 0) FillCalcHistograms(fHist->calc[4],calc,calc_par);
      else                     FillCalcHistograms(fHist->calc[5],calc,calc_par);

      if (fabs(calc_par->dtmin_crvc) < 30) {
        FillCalcHistograms(fHist->calc[6],calc,calc_par);
        if (calc->DiskID() == 0) FillCalcHistograms(fHist->calc[7],calc,calc_par);
        else                     FillCalcHistograms(fHist->calc[8],calc,calc_par);

        if (calc->Energy() > 20) {
          n_good_tc += 1;
        }
      }
    }
  }
//-----------------------------------------------------------------------------
// fill track histograms
//-----------------------------------------------------------------------------
  fHist->h_ntrk[0]->Fill(fNTrk);
  if (n_good_tc > 0) {
    fHist->h_ntrk[1]->Fill(fNTrk);
  }

  for (int i=0; i<fNTrk; i++) {
    TStrTrack* trk = fTrackBlock->Track(i);
    trk_param_t* tp = &fListOfTrkParam[i];
    FillTrkHistograms(fHist->trk[0],trk,tp);
    if (tp->intime) {
      
      FillTrkHistograms(fHist->trk[1],trk,tp);
      if (tp->calc->DiskID() == 0) {
        GetHeaderBlock()->Print("in-time disk0");
        FillTrkHistograms(fHist->trk[2],trk,tp);
      }
      else {
        GetHeaderBlock()->Print("in-time disk1");
        FillTrkHistograms(fHist->trk[3],trk,tp);
      }
    }
    
    if (tp->intime_calc) {
      FillTrkHistograms(fHist->trk[4],trk,tp);
      if (tp->calc->DiskID() == 0) {
        FillTrkHistograms(fHist->trk[5],trk,tp);
      }
      else if (tp->calc->DiskID() == 1) {
        FillTrkHistograms(fHist->trk[6],trk,tp);
      }
    }

    if ((trk->NHits() > 10) and (tp->intime_calc)) {
      FillTrkHistograms(fHist->trk[7],trk,tp);
      if (tp->calc->DiskID() == 0) {
        FillTrkHistograms(fHist->trk[8],trk,tp);
      }
      else if (tp->calc->DiskID() == 1) {
        FillTrkHistograms(fHist->trk[9],trk,tp);
      }
    }

    if ((trk->NHits() > 10) and (tp->intime_crvc)) {
      FillTrkHistograms(fHist->trk[10],trk,tp);
    }    
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::CalculateMissingTrkParameters() {
//-----------------------------------------------------------------------------
// extra track parameters
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNTrk; i1++) {
    TStrTrack*  trk = fTrackBlock->Track(i1);
    trk_param_t* tp = &fListOfTrkParam[i1];
//-----------------------------------------------------------------------------
// initialize parameter record
//-----------------------------------------------------------------------------
    tp->dtmin_tc    = 1.e6;
    tp->dtmin_crvc  = 1.e6;
    tp->dtmin_calc  = 1.e6;
    tp->tc          = nullptr;
    tp->calc        = nullptr;
    tp->crvc        = nullptr;
    tp->intime      = 0;
    tp->intime_calc = 0;
    tp->intime_crvc = 0;
    tp->dx_calc     = 1.e6;
    tp->dy_calc     = 1.e6;
    tp->xcrv        = 1.e6;
    tp->zcrv        = 1.e6;
    tp->dx_crvc     = 1.e6;
    tp->dz_crvc     = 1.e6;
//-----------------------------------------------------------------------------    
// find the closest time cluster - should always be there
//-----------------------------------------------------------------------------    
    for (int i2=0; i2<fNTc; i2++) {
      TStnTimeCluster* tc = fTcBlock->TimeCluster(i2);
      float dt = trk->fT0-tc->fT0;
      if (fabs(dt) < fabs(tp->dtmin_tc)) {
        tp->dtmin_tc = dt;
        tp->tc       = tc;
      }
    }
//-----------------------------------------------------------------------------
// determine the closest CRV coincidence
// as the calibration used time clusters, look at the time cluster T0
//-----------------------------------------------------------------------------
    float crv_time_offset = 0.; // 21; // today
    
    for (int i2=0; i2<fNCrvc; i2++) {
      TCrvCoincidenceCluster* crvc = fCrvcBlock->Cluster(i2);
      float dt = tp->tc->T0()-(crvc->StartTime()-crv_time_offset);
      if (fabs(dt) < fabs(tp->dtmin_crvc)) {
        tp->dtmin_crvc = dt;
        tp->crvc       = crvc;
      }
    }
    
    if (tp->crvc) {
      if (fabs(tp->dtmin_crvc) <  30) {
        tp->intime_crvc = 1;
      }
    
      float ycrv  = tp->crvc->Position()->Y() - 145.; // -80.;//-4280; // for kicks
      float dy    = ycrv-trk->fY0;
      
      tp->xcrv    = trk->fX0+(trk->fNx/trk->fNy)*dy;
      tp->zcrv    = trk->fZ0+(trk->fNz/trk->fNy)*dy;
      
      // X(CRVC) is defined in the global coordinate system
      float x_crvc = tp->crvc->Position()->X() + 3904.;
      tp->dx_crvc = tp->xcrv-x_crvc;
      
      tp->dz_crvc = tp->zcrv-(tp->crvc->Position()->Z()-23000.);  // approx
    }
//------------------------------;-----------------------------------------------
// determine the closest calorimeter cluster
//-----------------------------------------------------------------------------
    for (int i2=0; i2<fNCaloClusters; i2++) {
      TStnCluster* calc = fCaloClusterBlock->Cluster(i2);

      float dt = tp->tc->T0()-calc->Time();
      if (fabs(dt) < fabs(tp->dtmin_calc)) {
        tp->dtmin_calc = dt;
        tp->calc       = calc;
      }
    }
//-----------------------------------------------------------------------------
// extrapolate track to the closest cluster
//-----------------------------------------------------------------------------
    if (tp->calc) {
      double zc = 2383.; // 3560.; // guesswork for the disk0 position , not sure it is correct
      if (tp->calc->DiskID() == 1) {
        zc = 3517;
      }
      float dz    = zc-trk->fZ0;
      float tx    = trk->fX0+trk->fNx/trk->fNz*dz;
      float ty    = trk->fY0+trk->fNy/trk->fNz*dz;
      tp->dx_calc = tx-tp->calc->fX;
      tp->dy_calc = ty-tp->calc->fY;
    }

    if (fabs(tp->dtmin_calc) <  30) {
      tp->intime_calc = 1;
    }
    
    if ((fabs(tp->dtmin_tc  ) <  30) and
        (fabs(tp->dtmin_calc) <  30) and
        (fabs(tp->dtmin_crvc) <  30) and
        (tp->calc->Energy()   >= 20) and
        (trk->fNHits          >= 10)
        ) {
      tp->intime = 1;
    }
  }
  
  return 0;
}

//-----------------------------------------------------------------------------
int TDetTimeAnaModule::CalculateMissingParameters() {

  fNCalh10[0] = 0;
  fNCalh10[1] = 0;

  fListOfCalcParam.clear();
  if (fNCaloClusters > 0) fListOfCalcParam.resize(fNCaloClusters);
  
  fListOfTrkParam.clear();
  if (fNTrk > 0) fListOfTrkParam.resize(fNTrk);
  
  for (int i1=0; i1<fNCalh; i1++) {
    TCaloHit*  calh = fCaloHitBlock->Hit(i1);

    if (calh->fEDep > 10) {
      int disk = calh->Disk();
      fNCalh10[disk] += 1;
    }
  }
 
//-----------------------------------------------------------------------------
// extra parameters of the calorimeter clusters
//-----------------------------------------------------------------------------
  for (int i1=0; i1<fNCaloClusters; i1++) {
    TStnCluster*  calc = fCaloClusterBlock->Cluster(i1);

    calc_param_t& cp = fListOfCalcParam[i1];
//-----------------------------------------------------------------------------
// determine the closest time cluster
//-----------------------------------------------------------------------------
    cp.dtmin_tc = 1.e6;
    cp.tc       = nullptr;
    // determine the closest time cluster
    for (int i2=0; i2<fNTc; i2++) {
      TStnTimeCluster* tc = fTcBlock->TimeCluster(i2);
      float dt = calc->fTime-tc->fT0;
      if (fabs(dt) < fabs(cp.dtmin_tc)) {
        cp.dtmin_tc = dt;
        cp.tc       = tc;
      }
    }
//-----------------------------------------------------------------------------
// determine the closest CRV coincidence
//-----------------------------------------------------------------------------
    cp.dtmin_crvc = 1.e6;
    cp.crvc       = nullptr;

    float crv_time_offset = 0; // 21. // today
    
    for (int i2=0; i2<fNCrvc; i2++) {
      TCrvCoincidenceCluster* crvc = fCrvcBlock->Cluster(i2);
      float dt = calc->fTime-(crvc->StartTime()-crv_time_offset);
      if (fabs(dt) < fabs(cp.dtmin_crvc)) {
        cp.dtmin_crvc = dt;
        cp.crvc       = crvc;
      }
    }
//------------------------------;-----------------------------------------------
// determine the closest track
//-----------------------------------------------------------------------------
    cp.dtmin_trk = 1.e6;
    cp.trk       = nullptr;

    for (int i2=0; i2<fNTrk; i2++) {
      TStrTrack* trk = fTrackBlock->Track(i2);
      
      float dt = calc->fTime-trk->fT0;
      if (fabs(dt) < fabs(cp.dtmin_trk)) {
        cp.dtmin_trk = dt;
        cp.trk       = trk;
      }
    }
  }

  CalculateMissingTrkParameters();

  return 0;
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TDetTimeAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks
//-----------------------------------------------------------------------------
  RegisterDataBlock("CaloHitBlock"     ,"TCaloHitBlock"       ,&fCaloHitBlock     );
  RegisterDataBlock("CaloRecoDigiBlock","TCaloRecoDigiBlock"  ,&fCaloRecoDigiBlock);
  RegisterDataBlock("ClusterBlock"     ,"TStnClusterBlock"    ,&fCaloClusterBlock );
  RegisterDataBlock("CrvcBlock"        ,"TCrvClusterBlock"    ,&fCrvcBlock        );
  RegisterDataBlock("CrvpBlock"        ,"TCrvPulseBlock"      ,&fCrvpBlock        );
  RegisterDataBlock("TimeClusterBlock" ,"TStnTimeClusterBlock",&fTcBlock          );
  RegisterDataBlock("ComboHitBlock"    ,"TComboHitBlock"      ,&fChBlock          );
  RegisterDataBlock("TrackBlock"       ,"TStrTrackBlock"      ,&fTrackBlock       );
 
  return 0;
}


//_____________________________________________________________________________
int TDetTimeAnaModule::BeginRun() {
  fRunNumber = GetHeaderBlock()->RunNumber();
  TStntuple::Init(fRunNumber);
  fTpm->Init(fRunNumber);
//-----------------------------------------------------------------------------
// book histograms at begin run, there will be a flag telling whether we need
// to re-book them for every new run
// don't do that by default
//-----------------------------------------------------------------------------
  
  BookHistograms(fHist,fFolder);
  return 0;
}


//-----------------------------------------------------------------------------
// 2026-09-09
//-----------------------------------------------------------------------------
int TDetTimeAnaModule::Event(int IEntry) {

  //  TLorentzVector        mom;

  fTcBlock->GetEntry(IEntry);
  // don't need for the moment
  // fChBlock->GetEntry(IEntry);
  
  fCrvpBlock->GetEntry(IEntry);
  fCrvcBlock->GetEntry(IEntry);
  
  fCaloHitBlock->GetEntry(IEntry);
  fCaloRecoDigiBlock->GetEntry(IEntry);
  fCaloClusterBlock->GetEntry(IEntry);
  
  fTrackBlock->GetEntry(IEntry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNCrvc         = fCrvcBlock->NClusters();
  fNCrvp         = fCrvpBlock->NPulses();
  fNTc           = fTcBlock->NTimeClusters();
  fNCaloClusters = fCaloClusterBlock->NClusters();
  fNTrk          = fTrackBlock->NTracks();

  CalculateMissingParameters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TDetTimeAnaModule::Debug() {

  if (GetDebugBit(3) == 1) {
    // if (fNHelPos[1] > 0) {
    //   GetHeaderBlock()->Print(Form("N(CalHelixFinder helices hel > 0) = %2i",fNHelPos[1]));
    // }
  }

  if (GetDebugBit(4) == 1) {
    // if (fNHelNeg[1] > 0) {
    //   GetHeaderBlock()->Print(Form("N(CalHelixFinder helices hel < 0) = %2i",fNHelNeg[1]));
    // }
  }
}

//_____________________________________________________________________________
int TDetTimeAnaModule::EndJob() {
  return 0;
}


//-----------------------------------------------------------------------------
// do that for all ROCs and all FEBs
//-----------------------------------------------------------------------------
/*
# roc  feb   fit(crvp.time-tc.t0)  chi2dof
   1    1           530.724         3.004
   1    2           530.498         2.150
*/
int TDetTimeAnaModule::FitCaloTimeOffsets(float TMin, float TMax) {

  TH2F* h2 = fHist->h_dt_vs_sipmid;

  int nch = kNCaloChannels;
  for (int i=0; i<nch; i++) {
    
    std::string hpname = std::format("hpx_{:04d}",i);
    TH1D* hp = h2->ProjectionY(hpname.data(),i+1,i+1);
    fit_result_t* fr = &fFr[i];
    FitHistogram(hp,fr,TMin,TMax,10);
  }
  
  // done fitting, write results

  std::ofstream os("CaloTimeCalib_corr.txt");

  os << std::format("TABLE CalTimeCalib {}\n",fRunNumber);
  os << std::format("# corrections to  CalTimeCalib cid= XXX, to be subtracted \n");
  os << std::format("# sipmid     dT      sigma      chi2 \n");
  
  for (int i=0; i<nch; i++) {
    
    // to do no harm if no fit, dt should be initialized to zero 
    float dt{0}; 
    if (fFr[i].chi2dof > 0) {
      dt = fFr[i].p[1];
    }

    os << std::format("{:5}   {:8.3f}  {:8.3f}  {:8.3f}\n",i,dt,fFr[i].p[2],fFr[i].chi2dof);
  }
  os.close();
  
  return 0;
}

}
