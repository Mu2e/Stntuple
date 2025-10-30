///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TTriggerAnaModule_hh
#define Stntuple_ana_TTriggerAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TStnTriggerBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"
//-----------------------------------------------------------------------------
#include "Stntuple/ana/HistBase_t.h"
#include "Stntuple/ana/TrackPar_t.hh"
#include "Stntuple/ana/AnaDefs.hh"

namespace stntuple {
class TTriggerAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct TimeClusterHist_t : public stntuple::HistBase_t {
    TH1F* fTime;
    TH1F* fNHits;
    TH1F* fEnergy;
    TH1F* fR;
  };

  struct HelixHist_t : public stntuple::HistBase_t {
    TH1F* fNHits;
    TH1F* fNComboHits; 
    TH1F* fLambda;
    TH1F* fD0;
    TH1F* fNLoops;
    TH1F* fP;
    TH1F* fPt;
    TH1F* fT0;
    TH1F* fChi2XYNDof;
    TH1F* fChi2PhiZNDof;
    TH1F* fChi2TZNDof;
    TH1F* fSimpPDG1;
    TH1F* fMCP;
    TH1F* fDp;
    TH1F* fDpt;
    TH1F* fTZSlope;
    TH1F* fTZSlopeSig;
    TH1F* fHitRatio;
    TH1F* fPropDir;
  };

  struct TrackHist_t : public stntuple::HistBase_t {
    TH1F* fP;
    TH1F* fNActive;
    TH1F* fChi2Dof;
    TH1F* fT0;
    TH1F* fD0;
    TH1F* fZ0;
    TH1F* fTanDip;
    TH1F* fAlgMask;
  };

  struct TriggerHist_t : public stntuple::HistBase_t {
    TH1F* fBits;
    TH1F* fNPassedPaths;
    TH1F* fInstLum;
    TH1F* fNPassedPathInstLum[5];
    TH1F* fNPassedAprInstLum[3];
    TH1F* fNPassedCprInstLum;
    TH1F* fNPassedTprInstLum;
  };

  struct EventHist_t : public stntuple::HistBase_t {
    TH1F* fInstLum[2];
    TH1F* fRunNumber;
    TH1F* fEventNumber;
    TH1F* fNTimeClusters;
    TH1F* fNHelices   ;
    TH1F* fNTracks    ;
    TH1F* fMcMom      ;
    TH1F* fMcCosTh    ;
    TH1F* fTriggerBits;
    TH1F* fDtHelix;
    TH1F* fDpHelix;
    TH1F* fDptHelix;
    TH1F* fDradiusHelix;
    TH1F* fDd0Helix;
    TH1F* fDrcenterHelix;
    TH1F* fDtTrack;
    TH1F* fDpTrack;
    TH1F* fDptTrack;
    TH1F* fDz0Track;
    TH1F* fDd0Track;
    TH1F* fHelixPEff;
    TH1F* fHelixDeltaTEff;
    TH1F* fHelixD0Eff;
  };

//-----------------------------------------------------------------------------
  enum { kNTimeClusterHistSets = 100 };
  enum { kNHelixHistSets       = 100 };
  enum { kNTrackHistSets       = 200 };
  enum { kNTriggerHistSets     = 100 };
  enum { kNEventHistSets       = 100 };

  struct Hist_t {
    TimeClusterHist_t*  fTimeCluster [kNTimeClusterHistSets ];
    HelixHist_t*        fHelix    [kNHelixHistSets    ];
    TrackHist_t*        fTrack    [kNTrackHistSets    ];
    TriggerHist_t*      fTrigger  [kNTriggerHistSets  ];
    EventHist_t*        fEvent    [kNEventHistSets    ];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TGenpBlock*              fGenpBlock;
  TStnTimeClusterBlock*    fTimeClusterBlock;
  TStnHelixBlock*          fHelixBlock;
  TStnHelixBlock*          fAprHelixBlock;
  TStnTrackBlock*          fTrackBlock;
  TStnClusterBlock*        fClusterBlock;
  TStnTriggerBlock*        fTriggerBlock;

					// histograms filled
  Hist_t                   fHist;

  int                      fNGenp;
  int                      fNTracks;
  int                      fNTimeClusters;
  int                      fNHelices[2];
  int                      fNGoodTracks;
  stntuple::TrackPar_t     fTrackPar   [10];

  double                   fMinTrigMom;
  TGenParticle*            fParticle;
  int                      fPdgCode;
  int                      fProcessCode;
  double                   fMcMom;
  double                   fMcCosTh;
  double                   fWeight;

  std::vector<int>         fHelixIndexGood;
  std::vector<int>         fHelixIndexBad;
  std::vector<int>         fHelixIndexTail;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTriggerAnaModule(const char* name="Stntuple_TriggerAna", const char* title="Stntuple_TriggerAna");
  ~TTriggerAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void    SetPdgCode    (int Code   ) { fPdgCode     = Code ; }
  void    SetProcessCode(int Code   ) { fProcessCode = Code ; }
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookTimeClusterHistograms(stntuple::HistBase_t* Hist, const char* Folder);
  void    BookHelixHistograms      (stntuple::HistBase_t* Hist, const char* Folder);
  void    BookTrackHistograms      (stntuple::HistBase_t* Hist, const char* Folder);
  void    BookTriggerHistograms    (stntuple::HistBase_t* Hist, const char* Folder);
  void    BookEventHistograms      (stntuple::HistBase_t* Hist, const char* Folder);
  void    BookHistograms();

  void    FillTimeClusterHistograms(stntuple::HistBase_t* Hist, TStnTimeCluster*  TPeak);
  void    FillHelixHistograms      (stntuple::HistBase_t* Hist, TStnHelix*     Helix);
  void    FillTrackHistograms      (stntuple::HistBase_t* Hist, TStnTrack*     Trk);
  void    FillTriggerHistograms    (stntuple::HistBase_t* Hist);
  void    FillEventHistograms      (stntuple::HistBase_t* Hist, double Weight = 1);

  void    FillHistograms();

  void    InitTrackPar();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(stntuple::TTriggerAnaModule,0)
};
}
#endif
