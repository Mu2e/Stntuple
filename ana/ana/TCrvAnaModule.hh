///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_ana_TCrvAnaModule_hh__
#define __Stntuple_ana_TCrvAnaModule_hh__

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"
#include "Stntuple/obj/TCrvPulseBlock.hh"

namespace stntuple {
class TCrvAnaModule: public TStnModule {
public:

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fNCrvClusters;
    TH1F*    fNCrvCoincidences;
    TH1F*    fNCrvPulses;
  };

  struct GenpHist_t {
    TH1F*    fPdgCode[2];		// same distribution in different scale
    TH1F*    fGenID;			// 
    TH1F*    fZ0;			// 
    TH1F*    fT0;			// 
    TH1F*    fR0;			// 
    TH1F*    fP;			// 
    TH1F*    fCosTh;			// 
  };
					// histograms for the simulated CE
  struct SimpHist_t {
    TH1F*    fPdgCode;
    TH1F*    fMomTargetEnd;
    TH1F*    fMomTrackerFront;
  };

  struct CrvClusterHist_t {
    TH1F*    fSectorType;
    TH1F*    fNPulses;
    TH1F*    fNPe;
    TH1F*    fStartTime;
    TH1F*    fEndTime;
    TH1F*    fWidth;
    TH2F*    fXVsZ;
    TH2F*    fYVsZ;
  };

  struct CrvPulseHist_t {
    TH1F*    fNPe;
    TH1F*    fNPeHeight;
    TH1F*    fNDigis;
    TH1F*    fBar;
    TH1F*    fSipm;
    TH1F*    fTime;
    TH1F*    fHeight;
    TH1F*    fWidth;
    TH1F*    fChi2;
    TH1F*    fLeTime;
    TH1F*    fDt;
  };

  struct CrvCoincidenceHist_t {
    TH1F*    fSectorType;
    TH1F*    fNPulses;
  };

  enum { kNEventHistSets          = 100 };
  enum { kNGenpHistSets           = 100 };
  enum { kNSimpHistSets           = 100 };
  enum { kNCrvClusterHistSets     = 100 };
  enum { kNCrvPulseHistSets       = 100 };
  enum { kNCrvCoincidenceHistSets = 100 };

  struct Hist_t {
    EventHist_t*           fEvent         [kNEventHistSets];
    GenpHist_t*            fGenp          [kNGenpHistSets];
    SimpHist_t*            fSimp          [kNSimpHistSets];
    CrvClusterHist_t*      fCrvCluster    [kNCrvClusterHistSets];
    CrvCoincidenceHist_t*  fCrvCoincidence[kNCrvCoincidenceHistSets];
    CrvPulseHist_t*        fCrvPulse      [kNCrvPulseHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR

  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;
  TCrvClusterBlock*      fCrvClusterBlock;
  TCrvPulseBlock*        fCrvPulseBlock;

  // 					// additional track parameters (assume ntracks < 20)
  // TrackPar_t             fTrackPar[20];
					// histograms filled
  Hist_t                 fHist;
					// cut values
  double                 fPtMin;

  TSimParticle*          fSimp;		// electron or muon
  TGenParticle*          fParticle;		// electron or muon
  int                    fPdgCode;		// determines which one
  int                    fProcessCode;      

  int                    fNGenp;		// N(generated particles)

  double                 fMinT0;

  int                    fUseAllPulses;

  int                    fNCrvClusters;
  int                    fNCrvPulses;
  int                    fNCrvCoincidences;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TCrvAnaModule(const char* name="StntupleCrvAna", const char* title="Stntuple CrvAna");
  ~TCrvAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  void     SetProcessCode(int Code ) { fProcessCode = Code ; }
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
  void    BookGenpHistograms        (GenpHist_t*       Hist, const char* Folder);
  void    BookSimpHistograms        (SimpHist_t*       Hist, const char* Folder);
  void    BookCrvClusterHistograms  (CrvClusterHist_t* Hist, const char* Folder);
  void    BookCrvPulseHistograms    (CrvPulseHist_t*   Hist, const char* Folder);
  void    BookEventHistograms       (EventHist_t*      Hist, const char* Folder);

  void    FillEventHistograms       (EventHist_t*       Hist);
  void    FillGenpHistograms        (GenpHist_t*        Hist, TGenParticle* Genp );
  void    FillSimpHistograms        (SimpHist_t*        Hist, TSimParticle* Simp );
  void    FillCrvClusterHistograms  (CrvClusterHist_t*  Hist, TCrvCoincidenceCluster* CrvCl);
  void    FillCrvPulseHistograms    (CrvPulseHist_t*    Hist, TCrvRecoPulse* CrvCl);

  void    BookHistograms();
  void    FillHistograms();

  void    Debug();

  ClassDefOverride(stntuple::TCrvAnaModule,0)
};
}
#endif
