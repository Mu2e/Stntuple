///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_ana_TSpmcAnaModule_hh
#define Stntuple_ana_TSpmcAnaModule_hh

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "Stntuple/loop/TStnModule.hh"

#include "Stntuple/obj/TGenpBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"
#include "Stntuple/obj/TStepPointMCBlock.hh"

#include "Stntuple/base/TStnArrayI.hh"

#include "Stntuple/alg/TStntuple.hh"

#include "Stntuple/ana/HistBase_t.h"
#include "Stntuple/ana/SimPar_t.hh"
#include "Stntuple/ana/VDetData_t.hh"

#include "Stntuple/ana/AnaDefs.hh"

namespace stntuple {

class TSpmcAnaModule: public TStnModule {
public:
  enum { kMaxNSimp = 1000 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNSimp;
    TH1F*      fTMaxSimp[2];
    TH1F*      fTMaxSpmc;
    TH1F*      fLogWeight;
  };

  struct SimpHist_t : public HistBase_t {
    TH1F*      fVolumeID;		       //
    TH1F*      fStage;
    TH1F*      fGeneratorID;
    TH1F*      fTime;
    TH1F*      fParentPDG;
    TH1F*      fParentMom;
    TH1F*      fStartMom[2];

    TH2F*      fYVsX;
    TH2F*      fXEndVsZEnd;
    TH2F*      fYVsX_2480;
    TH2F*      fYVsX_2513;
    TH2F*      fCosThVsMom;
    TH2F*      fCosThVsMomPV;		// for antiprotons
  };

  struct StepPointMCHist_t : public HistBase_t {
    TH1F*      fVolumeID;		       //
    TH1F*      fGenIndex;		       //
    TH1F*      fSimID;
    TH1F*      fPDGCode[2];  // just different ranges
    TH1F*      fCreationCode;
    TH1F*      fParentSimID;
    TH1F*      fParentPDGCode;
    TH1F*      fEndProcessCode;

    TH1F*      fEDepTot;
    TH1F*      fEDepNio;
    TH1F*      fTime;			// in ns
    TH1F*      fTimeSec;		// in seconds, for R/A decays
    TH1F*      fStepLength;

    TH1F*      fMom[2];
    TH2F*      fCosThVsMom;
    TH1F*      fEKin;

    TH2F*      fYVsZ;
    TH2F*      fYVsX;
    TH2F*      fCosThVsMomPV;		// for antiprotons
  };

  struct VDetHist_t : public HistBase_t {
    TH1F*      fIndex   ;
    TH1F*      fPDGCode ;		       //
    TH1F*      fGenCode ;		       // generator code
    TH1F*      fMom[2]  ;
    TH1F*      fTime    ;
    TH1F*      fPTime   ;                // proper time
    TH2F*      fYVsX    ;                // different VD's have different orientation
    TH2F*      fYVsZ    ;                // fill both hist's
    TH1F*      fPt      ;                // transverse mom
    TH1F*      fPp      ;                // momentum component parallel to the solenoid axis
    TH1F*      fTanTh   ;		       // tan (pitch angle)
    TH1F*      fEKin    ;
    TH2F*      fCosThVsMom ;	       // cos (pitch angle) vs Mom
    TH2F*      fCosThVsMomPV;		// for antiprotons
    TH2F*      fTimeVsMom;
    TH2F*      fTimeVsMomW;
    TH2F*      fPTimeVsMom;
  };

  struct SimpData_t {
    int           fIndex;		// so far, not used
    TSimParticle* fParent;              // parent (in production vertex)
    int           fStage;               // for this particle
    double        fWeight;
  };

  struct SpmcData_t {
    TParticlePDG*  fParticle;		// so far, not used
  };

//-----------------------------------------------------------------------------
  enum { kNEventHistSets        =   100 };
  enum { kNStepPointMCHistSets  = 10000 };
  enum { kNSimpHistSets         = 10000 };
  enum { kNVDetHistSets         = 10000 };

  struct Hist_t {
    EventHist_t*        fEvent      [kNEventHistSets      ];
    SimpHist_t*         fSimp       [kNSimpHistSets       ];
    VDetHist_t*         fVDet       [kNVDetHistSets       ];
    StepPointMCHist_t*  fStepPointMC[kNStepPointMCHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TGenpBlock*           fGenpBlock;  
  TSimpBlock*           fSimpBlock;  
  TStepPointMCBlock*    fStepPointMCBlock;
  TStepPointMCBlock*    fVDetBlock;
  TStepPointMCBlock*    fSpmcBlockVDet;
					// histograms filled
  Hist_t                fHist;

  TString               fSpmcBlockName;
  TString               fVDetBlockName;

  TSimParticle*         fMuon;		// pointer to stopped muon (pend=0)
  TSimParticle*         fParent;
  TSimParticle*         fProton;

  TParticlePDG*         fParticleCache[5000];

  TDatabasePDG*         fPdgDb;

  int                   fNVDetHits  ;
  int                   fNVDet;
  stntuple::VDetData_t  fVDet[200];
  int                   fStageID;
  int                   fNSimp;
  int                   fStage;

  SimpData_t            fSimData[kMaxNSimp];

  TStntuple*            fStnt;
  double                fWeight;         // event weight, determined by the production cross section
  double                fTMaxSimp;	 // in seconds
  double                fTMaxSpmc;	 // in ns
					 // antiproton-specific : in the production vertex
  double                fPbarCosThPV;
  double                fPbarMomPV;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TSpmcAnaModule(const char* name="Stntuple_SpmcAna", const char* title="Stntuple SpmcAna");
  ~TSpmcAnaModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Hist_t*            GetHist        () { return &fHist;        }
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
  void SetSpmcBlockName(const char* Name) { fSpmcBlockName = Name; }
  void SetVDetBlockName(const char* Name) { fVDetBlockName = Name; }
  void SetStageID      (int ID) { fStageID = ID; }

  void          SetParticleCache(int PdgCode, TParticlePDG* P) { fParticleCache[2500+PdgCode] = P; }
  TParticlePDG* GetParticleCache(int PdgCode) { return fParticleCache[2500+PdgCode]; }
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
  void    BookEventHistograms        (HistBase_t* Hist, const char* Folder);
  void    BookSimpHistograms         (HistBase_t* Hist, const char* Folder);
  void    BookStepPointMCHistograms  (HistBase_t* Hist, const char* Folder);
  void    BookVDetHistograms         (HistBase_t* Hist, const char* Folder);

  void    FillEventHistograms        (HistBase_t* Hist);
  void    FillSimpHistograms         (HistBase_t* Hist, TSimParticle* Simp, SimpData_t* SimpData, double Weight = 1.);
  void    FillStepPointMCHistograms  (HistBase_t* Hist, TStepPointMC* Step, SpmcData_t* SpmcData, double Weight = 1.);
  void    FillVDetHistograms         (HistBase_t* Hist, TStepPointMC* Step,                       double Weight = 1.);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDefOverride(TSpmcAnaModule,0)
};

}
#endif
