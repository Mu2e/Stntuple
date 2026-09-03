///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TTrackSeedAnaModule_hh
#define Stntuple_ana_TTrackSeedAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

namespace stntuple {
class TTrackSeedAnaModule: public TStnModule {
public:
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1D*    fDioMom;
    TH1F*    fEleCosTh;
    TH1F*    fNHyp;
    TH1F*    fBestHyp[2];		// [0]: by chi2, [1]: by fit consistency
    TH1F*    fNGenp;                    // N(particles in GENP block)
    TH1F*    fNTrackSeeds;
  };

  struct TrackSeedHist_t {
    TH1F*    fNHits;	 
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fTanDip;   
    TH1F*    fChi2;
    TH1F*    fFitCons;
    TH1F*    fD0;
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
    TH1F*    fNStrawHits;
  };

  enum { kNEventHistSets   = 100 };
  enum { kNTrackSeedHistSets = 100 };
  enum { kNGenpHistSets    = 100 };
  enum { kNSimpHistSets    = 100 };

  struct Hist_t {
    EventHist_t*     fEvent     [kNEventHistSets];
    TrackSeedHist_t* fTrackSeed [kNTrackSeedHistSets];
    GenpHist_t*      fGenp      [kNGenpHistSets];
    SimpHist_t*      fSimp      [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used

  TStnTrackSeedBlock*    fTrackSeedBlock;
  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;
					// histograms filled
  Hist_t            fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fGeneratorCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNTrackSeeds;
  int               fNGenp;		// N(generated particles)


  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TTrackSeedAnaModule(const char* name="Stntuple_TrackSeedAna", const char* title="_Stntuple_ TrackSeedAna");
  ~TTrackSeedAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void               SetPdgCode      (int Code ) { fPdgCode       = Code ; }
  void               SetGeneratorCode(int Code ) { fGeneratorCode = Code ; }
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
  void    BookGenpHistograms    (GenpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms   (EventHist_t*   Hist, const char* Folder);
  void    BookSimpHistograms    (SimpHist_t*    Hist, const char* Folder);

  void    BookTrackSeedHistograms   (TrackSeedHist_t*   Hist, const char* Folder);

  void    FillEventHistograms    (EventHist_t* Hist);
  void    FillGenpHistograms     (GenpHist_t*    Hist, TGenParticle* Genp   );
  void    FillSimpHistograms     (SimpHist_t*    Hist, TSimParticle* Simp   );

  void    FillTrackSeedHistograms  (TrackSeedHist_t*   Hist, TStnTrackSeed*    TrkSeed);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDefOverride(TTrackSeedAnaModule,0)
};
}
#endif
