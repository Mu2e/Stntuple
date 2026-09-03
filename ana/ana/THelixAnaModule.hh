///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_THelixAnaModule_hh
#define Stntuple_ana_THelixAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TStnTimeClusterBlock.hh"
#include "Stntuple/obj/TStnHelixBlock.hh"
#include "Stntuple/obj/TStnTrackSeedBlock.hh"
#include "Stntuple/obj/TStnTrackBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"
#include "Stntuple/obj/TStrawHitBlock.hh"
#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/geom/TStnCrystal.hh"
#include "Stntuple/geom/TDiskCalorimeter.hh"

namespace  stntuple {
class THelixAnaModule: public TStnModule {
public:

  struct TrackPar_t {
    int     fNHPl;
    int     fNEPl;
    int     fNDPl;
    float   fDpF ;    // tracker-only resolution
    float   fDp0 ;
    float   fDp2 ;
    float   fDpFSt;
    double  fDioWt;

    double  fEcl;
    double  fEp;
    double  fDx;
    double  fDy;
    double  fDz;
    double  fDt;
    double  fDu;			// rotated residuals
    double  fDv;
    double  fChi2Match;
    double  fPath;
  };

  struct EventHist_t {
    TH1F*    fRv;			// MC truth information
    TH1F*    fZv;
    TH1F*    fEleMom;
    TH1F*    fEleCosTh;
    TH1F*    fNTimeClusters[2];
    TH1F*    fNHelices;
    TH1F*    fNHelPos;		       // N(helices) with helicity > 0
    TH1F*    fNHelNeg;		       // N(helices) with helicity < 0
    TH1F*    fNGenp;                   // N(particles in GENP block)
  };

  struct TimeClusterHist_t {
    TH1F*    fNHits;
    TH1F*    fNComboHits;
    TH1F*    fT0;
    TH1F*    fClusterEnergy;
  
  };


  struct HelixHist_t {
    TH1F*    fNHits;
    TH1F*    fHelicity;
    TH1F*    fClusterTime;
    TH1F*    fClusterEnergy;
    TH1F*    fRadius;    // fabs(1/omega)
    TH1F*    fMom;
    TH1F*    fPt;
    TH1F*    fLambda;   //dz/dphi
    TH1F*    fT0;
    TH1F*    fT0Err;
    TH1F*    fD0;
    TH1F*    fBestAlg;
    TH1F*    fAlgMask;
    TH1F*    fChi2XY;
    TH1F*    fChi2ZPhi;
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

  struct TrackEffHist_t {
    TH1F*    fPtMc;			// denominator
    TH1F*    fPtReco;			// numerator
  };

  enum { kNEventHistSets       = 100 };
  enum { kNHelixHistSets       = 300 };
  enum { kNTimeClusterHistSets = 200};
  enum { kNGenpHistSets        = 100 };
  enum { kNSimpHistSets        = 100 };

  struct Hist_t {
    EventHist_t*       fEvent       [kNEventHistSets];
    TimeClusterHist_t* fTimeCluster [kNTimeClusterHistSets];
    HelixHist_t*       fHelix       [kNHelixHistSets];
    GenpHist_t*        fGenp        [kNGenpHistSets];
    SimpHist_t*        fSimp        [kNSimpHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
					// 0: TPR, 1: CPR

  TStnTimeClusterBlock*  fTimeClusterBlock[2];
  TStnHelixBlock*        fHelixBlock;                 // 
  TStrawHitBlock*        fStrawHitBlock;
  TGenpBlock*            fGenpBlock;
  TSimpBlock*            fSimpBlock;

  TString                fHelixBlockName;

					// additional track parameters (assume ntracks < 20)
  TrackPar_t             fTrackPar[20];
					// histograms filled
  Hist_t                 fHist;
					// cut values
  double            fPtMin;

  TGenParticle*     fParticle;		// electron or muon
  int               fPdgCode;		// determines which one
  int               fProcessCode;      

  TSimParticle*     fSimp;
  double            fEleE;		// electron energy

  int               fNClusters;
  int               fNTimeClusters[2];

  int               fNHelices;	// 0: TPR, 1: CPR, 2: merged
  int               fNHelPos ;
  int               fNHelNeg ;

  int               fNStrawHits;
  int               fNGenp;		// N(generated particles)

  int               fNHyp;
  int               fBestHyp[10];
  int               fFillDioHist;

  double            fMinT0;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  THelixAnaModule(const char* name="Stntuple_HelixAna", const char* title="Stntuple_HelixAna");
  ~THelixAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*  GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  void     SetPdgCode    (int Code ) { fPdgCode     = Code ; }
  void     SetProcessCode(int Code ) { fProcessCode = Code ; }
  void     SetHelixBlockName(const char* Name) { fHelixBlockName = Name; }
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
  void    BookGenpHistograms        (GenpHist_t*    Hist, const char* Folder);
  void    BookSimpHistograms        (SimpHist_t*    Hist, const char* Folder);
  void    BookEventHistograms       (EventHist_t*   Hist, const char* Folder);

  void    BookTimeClusterHistograms (TimeClusterHist_t* Hist, const char* Folder);
  void    BookHelixHistograms       (HelixHist_t*       Hist, const char* Folder);

  void    FillEventHistograms       (EventHist_t* Hist);
  void    FillGenpHistograms        (GenpHist_t*  Hist, TGenParticle* Genp   );
  void    FillSimpHistograms        (SimpHist_t*  Hist, TSimParticle* Genp   );

  void    FillHelixHistograms       (HelixHist_t*       Hist, TStnHelix*       Helix      );
  void    FillTimeClusterHistograms (TimeClusterHist_t* Hist, TStnTimeCluster* TimeCluster);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDefOverride(THelixAnaModule,0)
};
}
#endif
