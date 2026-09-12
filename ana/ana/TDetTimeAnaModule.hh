///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_ana_TDetTimeAnaModule_hh__
#define __Stntuple_ana_TDetTimeAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TComboHitBlock.hh"

#include "Stntuple/obj/TCaloHitBlock.hh"
#include "Stntuple/obj/TCaloRecoDigiBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

#include "Stntuple/geom/TCrvChannelMap.hh"
#include "Stntuple/geom/TCaloChannelMap.hh"
#include "Stntuple/geom/TrkPanelMap.hh"

namespace stntuple {
class TDetTimeAnaModule: public TStnModule {
  
public:
  
  enum {
    kNStations         = 18,
    kNPlanes           = 36,
    kNPanelsPerStation = 12,
    kNCaloChannels     = 674*2*2,   // calorimeter
  };


  enum {
    kNEventHistSets      = 100,
    kNCrvdHistSets       = 100,
    kNCrvcHistSets       = 100,
    kNCrvpHistSets       = 100,
    kNRocHistSets        =  20,
    kNCalhHistSets       = 100,
    kNCalcHistSets       = 100,
    kNTrkHistSets        = 100,
  };
//-----------------------------------------------------------------------------
// indices
//-----------------------------------------------------------------------------
  struct TrkIndex_t {
    int sel;
    int slot;                 // 0-17
    int plane;                // offline
    int panel;                // offline
    int pnl12;                // panel index within the station (0-11)
    int mnid;
    int ch;
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
  
  struct CaloIndex_t {
    int sel   {-1};
    int crate {-1};
    int disk  {-1};
    int cid   {-1};                   //
  };
  
  struct calc_param_t {
    float                     dtmin_tc;           // from the closest TC
    TStnTimeCluster*          tc;
    float                     dtmin_crvc;     // from the closest CRVC
    TCrvCoincidenceCluster*   crvc;
    float                     dtmin_trk;                // from the closest track
    TStrTrack*                trk;                      // closest track
  };
                                                  // for now, a placeholder
  struct trk_param_t {
    int                       intime;
    int                       intime_calc;
    int                       intime_crvc;
    float                     dtmin_tc;           // from the closest TC
    TStnCluster*              calc;               // closest
    float                     dtmin_calc;         // from the closest CALC
    float                     dx_calc;
    float                     dy_calc;
    TStnTimeCluster*          tc;
    float                     dtmin_crvc;         // from the closest CRVC
    TCrvCoincidenceCluster*   crvc;
    float                     dxdz;
    float                     dydz;
    float                     xcrv;
    float                     zcrv;
    float                     dx_crvc;
    float                     dz_crvc;
  };
//-----------------------------------------------------------------------------
// histogram structures
//-----------------------------------------------------------------------------
  struct CalhHist_t {
    TH1F* h_dt;
  };
  
  struct ChannelHist_t {
    TH1F* h_ch;
    TH1F* h_dt;
  };

  struct CalcHist_t {
    TH1F*         h_edep;
    TH1F*         h_dt_tc;
    TH1F*         h_dt_crvc;
    TH2F*         h_dt_crvc_vs_dt_tc;
  };
  
  struct TrkHist_t {
    TH1F*         h_nhits;
    TH1F*         h_chi2d;
    TH1F*         h_t0;
    TH1F*         h_dt_crvc;
    TH1F*         h_dt_calc;
    TH1F*         h_dt_tc;
    TH1F*         h_dx_calc;
    TH1F*         h_dy_calc;
    TH2F*         h_dx_calc_vs_dxdz;
    TH2F*         h_dy_calc_vs_dydz;
    TH1F*         h_dxdz;
    TH1F*         h_dydz;
    TH1F*         h_xcrv;
    TH1F*         h_zcrv;
    TH1F*         h_dx_crvc;
    TH1F*         h_dz_crvc;
    TH2F*         h_dx_crvc_vs_dxdy;
    TH2F*         h_dz_crvc_vs_dzdy;
  };
  
  struct DiskHist_t {
    // BoardHist_t*  board[30];
    TH1F*         h_ch;
    TH2F*         h_board_vs_dt;
  };
  
  struct Hist_t {
    DiskHist_t* disk[2];                // 2 disks
    TH2F*       h_dt_vs_sipmid;
    TH1F*       h_sipmid;               // occupancy offline channel
    TH2F*       h_n2_vs_n1;
    TH1F*       h_ntrk[2];
    CalhHist_t* calh   [kNCalhHistSets];
    CalcHist_t* calc   [kNCalcHistSets];
    TrkHist_t*  trk    [kNTrkHistSets ];
  };


//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR
  TStnTimeClusterBlock*     fTcBlock;
  TCrvClusterBlock*         fCrvcBlock;
  TCrvPulseBlock*           fCrvpBlock;
  TComboHitBlock*           fChBlock;
  TStrTrackBlock*           fTrackBlock;
  TStnClusterBlock*         fCaloClusterBlock;
  TCaloHitBlock*            fCaloHitBlock;
  TCaloRecoDigiBlock*       fCaloRecoDigiBlock;

  TrkPanelMap*              fTpm;
  TCaloChannelMap*          fCaloChannelMap;
  TCrvChannelMap*           fCrvChannelMap;
  
  int                       fRunNumber;

  int                       fMaxEvent;   // for X-axis truncation
  int                       fNEvents;
 					
  Hist_t*                   fHist;       // histograms to be filled

  int                       fNCalh10[2]; // N(hits above 10 MeV)
  int                       fNTrk;
  int                       fNCalh;
  int                       fNCalrd;
  int                       fNCaloClusters;
  int                       fNTc;
  int                       fNCrvc;
  int                       fNCrvp;
  int                       fNCrvd;

  std::vector<calc_param_t> fListOfCalcParam;
  std::vector<trk_param_t>  fListOfTrkParam;

  fit_result_t              fFr[kNCaloChannels];
  fit_result_t*             fFrRef;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDetTimeAnaModule(const char* name="DetTimeAna", const char* title="Stntuple DetTimeAna");
  ~TDetTimeAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  // void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  // void     SetProcessCode(int Code ) { fProcessCode = Code ; }
  //  void     SetUseAllPulses(int Flag) { fUseAllPulses = Flag; }
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
  int              BookCalcHistograms (CalcHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookCalhHistograms (CalhHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookDiskHistograms (DiskHist_t*   Hist, CaloIndex_t* Index, TFolder* Folder);
  int              BookTrkHistograms  (TrkHist_t*    Hist, TrkIndex_t*  Index, TFolder* Folder);
  int              BookHistograms     (Hist_t*       Hist, TFolder*     Folder);

  int              CalculateMissingParameters   ();
  int              CalculateMissingTrkParameters();
  
  int              FillCalcHistograms (CalcHist_t* Hist, TStnCluster*   Calc, calc_param_t* Cp);
  int              FillDiskHistograms (DiskHist_t* Hist, TCaloRecoDigi* Calrd);
  int              FillTrkHistograms  (TrkHist_t*  Hist, TStrTrack*     Trk , trk_param_t* Tp);
  int              FillHistograms     ();

  int              FitCaloTimeOffsets (float TMin = 1., float TMax = -1.);
  //  int              PrintTimeCorrections();

  void             Debug();

  ClassDefOverride(stntuple::TDetTimeAnaModule,0)
};
}
#endif
