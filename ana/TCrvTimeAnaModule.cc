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
#include "ana/TCrvTimeAnaModule.hh"


ClassImp(stntuple::TCrvTimeAnaModule)

namespace stntuple {
//-----------------------------------------------------------------------------
TCrvTimeAnaModule::TCrvTimeAnaModule(const char* name, const char* title):
  TStnModule(name,title)
{
  // fPtMin = 1.;
  // fMinT0 = 0; // do not cut on time by default
  fHist = new Hist_t;
}

//-----------------------------------------------------------------------------
TCrvTimeAnaModule::~TCrvTimeAnaModule() {
}

//-----------------------------------------------------------------------------
  void TCrvTimeAnaModule::BookEventHistograms(EventHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {
  
  std::string prefix = std::format("run:{:06d}",fRunNumber);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  name  = "feb_vs_ch";
  title = std::format("{} : FEB vs CH",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

  name  = "dt_vs_sbid";
  title = std::format("{} : dt vs SBID",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_sbid,name.data(),title.data(),1500,0,1500,200,0,1000,Folder);

  name  = "feb_vs_sbid_0";
  title = std::format("{} : dt vs SBID ROC=1",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_sbid[0],name.data(),title.data(),600,0,600,30,0,30,Folder);

  name  = "feb_vs_sbid_1";
  title = std::format("{} : dt vs SBID ROC=2",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_sbid[1],name.data(),title.data(),600,0,600,30,0,30,Folder);

  name  = "och";
  title = std::format("{} : offline channel ID",prefix);
  fBookHist->HBook1F(Hist->h_och,name.data(),title.data(),2400,0,2400,Folder);
}
  
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::BookCrvcHistograms(CrvcHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{:02d}",fRunNumber,Index->sel);
  std::string name, title;

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBookHist->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...
}

 
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::BookCrvpHistograms(CrvpHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  std::string prefix = std::format("run:{:06d} sel:{} roc:{:02d} feb:{}",
                                   fRunNumber,Index->sel, Index->roc, Index->feb);
  std::string name, title;

  name  = "ph";
  title = std::format("{} : ph",prefix);
  fBookHist->HBook1F(Hist->h_ph,name.data(),title.data(),100,0,1000,Folder);   // in us...

  name  = "npes";
  title = std::format("{} : npes",prefix);
  fBookHist->HBook1F(Hist->h_npes,name.data(),title.data(),100,0,500,Folder);   // in us...

  name  = "time";
  title = std::format("{} : time",prefix);
  fBookHist->HBook1F(Hist->h_time,name.data(),title.data(),100,0,1.e5,Folder);   // in us...

  name  = "dt";
  title = std::format("{} : dt",prefix);
  fBookHist->HBook1F(Hist->h_dt,name.data(),title.data(),1000,-1000,1000,Folder);   // in us...

  name  = "feb";
  title = std::format("{} : feb",prefix);
  fBookHist->HBook1F(Hist->h_feb,name.data(),title.data(),100,0,100,Folder);   // in us...

  name  = "ch";
  title = std::format("{} : ch",prefix);
  fBookHist->HBook1F(Hist->h_ch,name.data(),title.data(),2000,0,2000,Folder);   // in us...

  name  = "dt_vs_feb";
  title = std::format("{} : dt vs feb",prefix);
  fBookHist->HBook2F(Hist->h_dt_vs_feb,name.data(),title.data(),100,0,100,1000,-1000,1000,Folder);   // in us...
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::BookFebHistograms(FebHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  // Index_t index;

  std::string prefix = std::format("run:{:06d} roc:{} feb:{:02d}",fRunNumber,Index->roc, Index->feb);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  // name  = "dt";
  // title = std::format("{} : T(pulse)-T(trk TC)",prefix);
  // fBookHist->HBook1F(Hist->h_dt,name.data(),title.data(),400,-1000,1000,Folder);

  name  = "ch_vs_dt";
  title = std::format("{} : channel vs [T(pulse)-T(trk TC)]",prefix);
  fBookHist->HBook2F(Hist->h_ch_vs_dt,name.data(),title.data(),400,-1000,1000,64,0,64,Folder);

  // name  = "feb_vs_ch";
  // title = std::format("{} : FEB vs CH",prefix);
  // fBookHist->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);
}



//-----------------------------------------------------------------------------
// by FEB: 32 scintillation bars per FEB 
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::BookRocHistograms(RocHist_t* Hist, CrvIndex_t* Index, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  // Index_t index;

  std::string prefix = std::format("run:{:06d} roc:{:02d}",fRunNumber,Index->roc);
  std::string name, title;

  name  = "sbid";
  title = std::format("{} : SBID",prefix);
  fBookHist->HBook1F(Hist->h_sbid,name.data(),title.data(),1500,0,1500,Folder);

  name  = "feb_vs_dt";
  title = std::format("{} : FEB vs_dt",prefix);
  fBookHist->HBook2F(Hist->h_feb_vs_dt,name.data(),title.data(),1000,-1000,1000,30,0,30,Folder);

  // name  = "feb_vs_ch";
  // title = std::format("{} : FEB vs CH",prefix);
  // fBookHist->HBook2F(Hist->h_feb_vs_ch,name.data(),title.data(),64,0,64,30,0,30,Folder);

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  std::string folder_name = std::format("crvp");
  TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
  if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
  Hist->crvp = new CrvpHist_t();
  BookCrvpHistograms(Hist->crvp,Index,fol);

//-----------------------------------------------------------------------------
// by FEB - 24 or 25 FEBs per ROC ?
//-----------------------------------------------------------------------------
  int n_feb_histsets(30);

  for (int i=0; i<n_feb_histsets; i++) {
    std::string folder_name = std::format("feb_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->feb[i] = new FebHist_t();
    Index->feb = i;
    BookFebHistograms(Hist->feb[i],Index,fol);
  }
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::BookHistograms(Hist_t* Hist, TFolder* Folder) {

  // std::string prefix = std::format("");
  // std::string name, title;

  CrvIndex_t index;
//-----------------------------------------------------------------------------
// Event histograms
//-----------------------------------------------------------------------------
  int book_event_histset[kNEventHistSets];

  for (int i=0; i<kNEventHistSets; i++) { book_event_histset[i] = 0; }

  book_event_histset[0] = 1;

  for (int i=0; i<kNEventHistSets; i++) {
    if (book_event_histset[i] == 0) continue;
    std::string folder_name = std::format("evt_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->event[i] = new EventHist_t();
    BookEventHistograms(Hist->event[i],&index,fol);
  }


  BookEventHistograms(Hist->event[0],&index,Folder);
//-----------------------------------------------------------------------------
// CRV coincidence clusters
//-----------------------------------------------------------------------------
  int book_crvc_histset[10];
  int n_crvc_histsets(10);

  for (int i=0; i<n_crvc_histsets; i++) { book_crvc_histset[i] = 0; }

  book_crvc_histset[0] = 1;

  for (int i=0; i<n_crvc_histsets; i++) {
    if (book_crvc_histset[i] == 0) continue;
    std::string folder_name = std::format("crvc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvc[i] = new CrvcHist_t();
    BookCrvcHistograms(Hist->crvc[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// CRV reco pulses
//-----------------------------------------------------------------------------
  int book_crvp_histset[10];
  int n_crvp_histsets(10);

  for (int i=0; i<n_crvp_histsets; i++) { book_crvp_histset[i] = 0; }

  book_crvp_histset[0] = 1;             // all
  book_crvp_histset[1] = 1;             // ntc=1, dt_530 < 30
  book_crvp_histset[2] = 1;             // ntc=1, dt_590 < 30

  for (int i=0; i<n_crvp_histsets; i++) {
    if (book_crvp_histset[i] == 0) continue;
    std::string folder_name = std::format("crvp_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->crvp[i] = new CrvpHist_t();
    index.sel = i;
    BookCrvpHistograms(Hist->crvp[i],&index,fol);
  }

//-----------------------------------------------------------------------------
// by ROC, links 0 and 3 --> rocs #1 and #4
//-----------------------------------------------------------------------------
  int book_roc_histset[10];
  int n_roc_histsets(10);

  for (int i=0; i<n_roc_histsets; i++) { book_roc_histset[i] = 0; }

  book_roc_histset[1] = 1;
  book_roc_histset[2] = 1;
  book_roc_histset[4] = 1;

  for (int i=0; i<n_roc_histsets; i++) {
    if (book_roc_histset[i] == 0) continue;
    std::string folder_name = std::format("roc_{:02d}",i);
    TFolder* fol = (TFolder*) Folder->FindObject(folder_name.data());
    if (! fol) fol = Folder->AddFolder(folder_name.data(),folder_name.data());
    Hist->roc[i] = new RocHist_t();
    index.roc = i;
    BookRocHistograms(Hist->roc[i],&index,fol);
  }
}


//-----------------------------------------------------------------------------
// need MC truth branch
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::FillEventHistograms(EventHist_t* Hist) {
  //  double            cos_th(-2.), p(-1.);
  // double            xv(-1.e6), yv(-1.e6), rv(-1.e6), zv(-1.e6);
  // TLorentzVector    mom;

  // if (fParticle) {
  //   fParticle->Momentum(mom);
  //   //    p      = mom.P();
  //   //    cos_th = mom.Pz()/p;
  //   xv = fParticle->Vx()+3904.;
  //   yv = fParticle->Vy();
  //   rv = sqrt(xv*xv+yv*yv);
  //   zv = fParticle->Vz();
  // }

  // //  Hist->fEleMom->Fill(p);
  // //  Hist->fEleCosTh->Fill(cos_th);
  // Hist->fRv->Fill(rv);
  // Hist->fZv->Fill(zv);

  // Hist->fNCrvClusters->Fill(fNCrvClusters);
  // Hist->fNCrvCoincidences->Fill(fNCrvCoincidences);
  // Hist->fNCrvPulses->Fill(fNCrvPulses);
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::FillCrvcHistograms(CrvcHist_t* Hist, TCrvCoincidenceCluster* CrvCluster) {

  // Hist->fSectorType->Fill(CrvCluster->SectorType());
  // Hist->fNPulses->Fill(CrvCluster->NPulses());
  // Hist->fNPe->Fill(CrvCluster->NPe());
  // Hist->fStartTime->Fill(CrvCluster->StartTime());
  // Hist->fEndTime->Fill(CrvCluster->EndTime());

  // float width  = CrvCluster->EndTime()-CrvCluster->StartTime();
  // Hist->fWidth->Fill(width);

  // float x = CrvCluster->Position()->X();
  // float y = CrvCluster->Position()->Y();
  // float z = CrvCluster->Position()->Z();

  // Hist->fXVsZ->Fill(z,x);
  // Hist->fYVsZ->Fill(z,y);
}
//-----------------------------------------------------------------------------
// need to optimize the filling time
//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::FillCrvpHistograms(CrvpHist_t* Hist, TCrvRecoPulse* Crvp) {
  // filling histograms: plot time differences between
  Hist->h_ph->Fill(Crvp->fPh);
  Hist->h_npes->Fill(Crvp->Pes());
  Hist->h_time->Fill(Crvp->Time());
  Hist->h_feb->Fill(Crvp->Feb());
  Hist->h_ch->Fill(Crvp->FebCh());
}

//-----------------------------------------------------------------------------
// register data blocks and book histograms
//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::BeginJob() {
//-----------------------------------------------------------------------------
// register data blocks, 'HelixBlock' - OR of TPR and CPR
//-----------------------------------------------------------------------------
  RegisterDataBlock("CrvcBlock"       ,"TCrvClusterBlock"    ,&fCrvcBlock);
  RegisterDataBlock("CrvpBlock"       ,"TCrvPulseBlock"      ,&fCrvpBlock);
  RegisterDataBlock("TimeClusterBlock","TStnTimeClusterBlock",&fTcBlock  );
  RegisterDataBlock("ComboHitBlock"   ,"TComboHitBlock"      ,&fChBlock  );
//-----------------------------------------------------------------------------
// book histograms
//-----------------------------------------------------------------------------
  BookHistograms(fHist,fFolder);

  return 0;
}


//_____________________________________________________________________________
int TCrvTimeAnaModule::BeginRun() {
  int rn = GetHeaderBlock()->RunNumber();
  TStntuple::Init(rn);
  return 0;
}


//_____________________________________________________________________________
void TCrvTimeAnaModule::FillHistograms() {

  EventHist_t* ehr = fHist->event[0];
  
  for (int i2=0; i2<fNCrvp; i2++) {
    TCrvRecoPulse*  crvp = fCrvpBlock->Pulse(i2);

    // this is a global histogram
    ehr->h_feb_vs_ch->Fill(crvp->fFebCh,crvp->fFeb);
                                        // have two histograms - one per roc , to color them
                                        // ROCs 1 and 2 --> hists 0 and 1
    ehr->h_feb_vs_sbid[crvp->fRoc-1]->Fill(crvp->fSbid,crvp->fFeb);
    ehr->h_sbid->Fill(crvp->fSbid);

    RocHist_t* rhr = fHist->roc[crvp->fRoc];
    rhr->h_sbid->Fill(crvp->fSbid);

    FebHist_t* fhr = rhr->feb[crvp->fFeb];
    fhr->h_sbid->Fill(crvp->fSbid);
                                        // this is a global histogram
    int och = crvp->OfflineChID();
    ehr->h_och->Fill(och);
//-----------------------------------------------------------------------------
// CRVP[0] : all pulses
//-----------------------------------------------------------------------------
    FillCrvpHistograms(fHist->crvp[0],crvp);
  }
//-----------------------------------------------------------------------------
// double-nested loops start here
//-----------------------------------------------------------------------------
  for (int i=0; i<fNTc; i++) {
    TStnTimeCluster*  tc = fTcBlock->TimeCluster(i);

    for (int i2=0; i2<fNCrvc; i2++) {
      TCrvCoincidenceCluster*  crvc = fCrvcBlock->Cluster(i2);
      float dt = crvc->StartTime()-tc->T0();
      
      fHist->crvc[0]->h_dt->Fill(dt);
    }

    for (int i2=0; i2<fNCrvp; i2++) {
      TCrvRecoPulse*  crvp = fCrvpBlock->Pulse(i2);
      float dt       = crvp->Time()-tc->T0();
      //int feb = crvp->feb;

      CrvpHist_t* crvp_hr = fHist->crvp[0];
      
      crvp_hr->h_dt->Fill(dt);
      ehr->h_dt_vs_sbid->Fill(crvp->Sbid(),dt);

      RocHist_t* roc_hr = fHist->roc[crvp->fRoc];
      roc_hr->h_feb_vs_dt->Fill(dt,crvp->fFeb);

      FebHist_t* feb_hr = fHist->roc[crvp->fRoc]->feb[crvp->fFeb];
      feb_hr->h_ch_vs_dt->Fill(dt,crvp->fFebCh);

      if (fNTc == 1) {
        if      (fabs(dt - 530) < 30) {
//-----------------------------------------------------------------------------
// CRVP[1] : first peak
//-----------------------------------------------------------------------------
          FillCrvpHistograms(fHist->crvp[1],crvp);
        }
        else if (fabs(dt - 590) < 30) {
//-----------------------------------------------------------------------------
// CRVP[2] : second peak
//-----------------------------------------------------------------------------
          FillCrvpHistograms(fHist->crvp[2],crvp);
        }
      }
    }
  }
}



//-----------------------------------------------------------------------------
// 2014-04-30: it looks that reading the straw hits takes a lot of time - 
//              turn off by default by commenting it out
//-----------------------------------------------------------------------------
int TCrvTimeAnaModule::Event(int ientry) {

  //  TLorentzVector        mom;

  fTcBlock->GetEntry(ientry);
  fChBlock->GetEntry(ientry);
  fCrvpBlock->GetEntry(ientry);
  fCrvcBlock->GetEntry(ientry);
//-----------------------------------------------------------------------------
// assume electron in the first particle, otherwise the logic will need to 
// be changed
//-----------------------------------------------------------------------------
  fNCrvc = fCrvcBlock->NClusters();
  fNCrvp = fCrvpBlock->NPulses();
  fNTc   = fTcBlock->NTimeClusters();
  
  FillHistograms();

  Debug();

  return 0;		       
}

//-----------------------------------------------------------------------------
void TCrvTimeAnaModule::Debug() {

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
int TCrvTimeAnaModule::EndJob() {
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
int TCrvTimeAnaModule::FitFebTimeOffsets(float TMin, float TMax) {

  for (int i=1; i<3; i++) {
    TH2F* h2 = fHist->roc[i]->h_feb_vs_dt;
    
    // fit Y-slices, use Ralf's integers
    for (int j=1; j<25; j++) {
      std::string hpname = std::format("hpx_{:02d}",j);
      TH1D* hp = h2->ProjectionX(hpname.data(),j+1,j+1);
      fit_result_t* fr = &fFr[i][j];
      FitHistogram(hp,fr,TMin,TMax,100);
    }
  }
  
  // done fitting, print results

  std::cout << std::format("# roc  feb   fit(crvp.time-tc.t0) dT(i-0)    sigma_i       chi2dof\n");
  
// find the first converged fit

  fFrRef = nullptr;

  bool ref_ch_found(false);
  
  for (int i=1; i<3; i++) {
    for (int j=1; j<25; j++) {
      if (fFr[i][j].chi2dof > 0) {
        fFrRef = &fFr[i][j];
        std::cout << std::format("reference channel: roc:{:2} feb:{:2} dt0:{:8.3f}\n",i,j,fFrRef->p[1]);
        ref_ch_found = true;
        break;
      }
    }
    if (ref_ch_found) break;
  }

  for (int i=1; i<3; i++) {
    for (int j=1; j<25; j++) {
      fit_result_t* fr = &fFr[i][j] ;
      float dt(0);
      if (fr->chi2dof > 0) {
        dt = fr->p[1]-fFrRef->p[1];
      }
      std::cout << std::format(" {:2d} {:4d}      {:9.3f}      {:9.3f}   {:9.3f}    {:9.3f}\n",
                               i,j,fr->p[1],dt,fr->p[2],fr->chi2dof);
    }
  }
  return 0;
}

//-----------------------------------------------------------------------------
// corrections are aimed to align in time all FEBs with FEB[1][1]
// to be called AFTER FitFebTimeOffset - that defines fFrRef
//------------------------------------------------------------------------------
int TCrvTimeAnaModule::PrintTimeCorrections() {
  
  std::ofstream os("CrvTime_corr.txt");

  float dt0 = fFrRef->p[1];  // roc=1 feb=1

  for (int i=0; i<2304; i++) {
    
    TCrvChannelMap::Data_t* dat = fCcm->ch_data_by_offline(i);

    float dt{0};
    if (dat != nullptr) {
      if (fFr[dat->roc][dat->feb].chi2dof > 0) {
        // dont correct FEBs with no fit
        dt = fFr[dat->roc][dat->feb].p[1] - dt0; // should be initialized to zero
      }
    }

    os << std::format("{:5}   {:8.3f}\n",i,dt);
  }

  os.close();
  
  return 0;
}

}
