///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_ana_TCrvTimeAnaModule_hh__
#define __Stntuple_ana_TCrvTimeAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

// #include "Stntuple/obj/TGenpBlock.hh"
// #include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TComboHitBlock.hh"

#include "Stntuple/geom/TCrvChannelMap.hh"

namespace stntuple {
class TCrvTimeAnaModule: public TStnModule {
  
public:

  enum {
    kNEventHistSets      = 100,
    kNCrvdHistSets       = 100,
    kNCrvcHistSets       = 100,
    kNCrvpHistSets       = 100,
    kNRocHistSets        =  20,
  };

  struct CrvIndex_t {
    int sel  {-1};
    int sbid {-1};
    int sipm {-1};
    int och  {-1};                           // offline channel - 4*sbid+sipm
    int roc  {-1};                           // 1-18 ??? 
    int feb  {-1};                           // 1-24 in offline domain
    int ch   {-1};                           // channel within the FEB (0-63)
  };
  
  struct CrvcHist_t {
    TH1F* h_dt;
  };
  
  struct ChannelHist_t {
    TH1F* h_ch;
    TH1F* h_dt;
  };

  struct CrvpHist_t {
    TH1F* h_ph;
    TH1F* h_npes;
    TH1F* h_time;
    TH1F* h_dt;
    TH1F* h_feb;
    TH1F* h_ch;
    TH2F* h_dt_vs_feb;
  };

  struct FebHist_t {
    ChannelHist_t ch[64];
    CrvpHist_t*   crvp;               // for the whole FEB
    TH1F*         h_sbid;
    TH2F*         h_ch_vs_dt;
  };
  
  struct RocHist_t {
    CrvpHist_t*   crvp;               // for the whole ROC
    FebHist_t*    feb [30];
    TH1F*         h_sbid;
    TH2F*         h_feb_vs_dt;
  };

  struct EventHist_t {
    TH1F*       h_sbid;
    TH2F*       h_feb_vs_ch;
    TH2F*       h_dt_vs_sbid;
    TH2F*       h_feb_vs_sbid[2];       // one per ROC
    TH1F*       h_och;                  // occupancy offline channel
  };
  
  struct Hist_t {
    EventHist_t* event[kNEventHistSets];
    // CrvdHist_t*  crvd [kNCrvdHistSets];
    CrvcHist_t*  crvc [kNCrvcHistSets];
    CrvpHist_t*  crvp [kNCrvpHistSets];
    RocHist_t*   roc  [kNRocHistSets];
  };

//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR
  TStnTimeClusterBlock* fTcBlock;
  TCrvClusterBlock*     fCrvcBlock;
  TCrvPulseBlock*       fCrvpBlock;
  TComboHitBlock*       fChBlock;

  // TrkPanelMap_t*    fTpm;
  
  TCrvChannelMap*   fCcm;
  int               fRunNumber;

  int               fMaxEvent;             // for X-axis truncation
  int               fNEvents;
 					// histograms filled
  Hist_t*           fHist;
		           	// cut values
  // double            fMinT0;

  int               fUseAllPulses;

  int               fNCrvc;
  int               fNCrvp;
  int               fNTc;               // N(time clusters)

  fit_result_t      fFr[10][25];
  fit_result_t*     fFrRef;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCrvTimeAnaModule(const char* name="CrvTimeAna", const char* title="Stntuple CrvTimeAna");
  ~TCrvTimeAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  // void     SetProcessCode(int Code ) { fProcessCode = Code ; }
  void     SetUseAllPulses(int Flag) { fUseAllPulses = Flag; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int     BeginJob() override;
  virtual int     BeginRun() override;
  virtual int     Event   (int ientry) override;
  virtual int     EndJob  () override;
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookCrvcHistograms  (CrvcHist_t*   Hist, CrvIndex_t* Index, TFolder* Folder);
  void    BookCrvpHistograms  (CrvpHist_t*   Hist, CrvIndex_t* Index, TFolder* Folder);
  void    BookEventHistograms (EventHist_t*  Hist, CrvIndex_t* Index, TFolder* Folder);
  void    BookFebHistograms   (FebHist_t*    Hist, CrvIndex_t* Index, TFolder* Folder);
  void    BookRocHistograms   (RocHist_t*    Hist, CrvIndex_t* Index, TFolder* Folder);
  void    BookHistograms      (Hist_t* Hist, TFolder* Folder);

  void    FillEventHistograms (EventHist_t*  Hist);
  void    FillCrvcHistograms  (CrvcHist_t*   Hist, TCrvCoincidenceCluster* Crvc);
  void    FillCrvpHistograms  (CrvpHist_t*   Hist, TCrvRecoPulse*          Crvp);
  void    FillHistograms();

  
  // if TMax > TMin, use them as limits, otherwise determine them automatically
  int     FitFebTimeOffsets   (float TMin = 1., float TMax = -1.);
  
  int     PrintTimeCorrections();
  

  void    Debug();

  ClassDefOverride(stntuple::TCrvTimeAnaModule,0)
};
}
#endif
