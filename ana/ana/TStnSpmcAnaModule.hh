///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef murat_ana_TStnSpmcAnaModule_hh
#define murat_ana_TStnSpmcAnaModule_hh

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

// #include "murat/ana/HistBase_t.h"
// #include "murat/ana/SimPar_t.hh"
#include "Stntuple/ana/StnVDetData_t.hh"

class TStnSpmcAnaModule: public TStnModule {
public:
  enum { kMaxNSimp = 1000 };
//-----------------------------------------------------------------------------
//  histograms
//-----------------------------------------------------------------------------
  struct EventHist_t { // : public HistBase_t {
    TH1F*      fRunNumber;
    TH1F*      fEventNumber;
    TH1F*      fNSimp;
    TH1F*      fNSpmc;
    TH1F*      fNSpmc150;
  };

  struct SimpHist_t { // : public HistBase_t {
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

  struct SpmcHist_t { // : public HistBase_t {
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
    TH1F*      fTime;
    TH1F*      fStepLength;

    TH1F*      fMom[2];
    TH2F*      fCosThVsMom;
    TH1F*      fEKin;

    TH2F*      fYVsZ;
    TH2F*      fYVsX;
    TH2F*      fCosThVsMomPV;		// for antiprotons
  };

  struct VDetHist_t { // : public HistBase_t {
    TH1F*      fIndex   ;
    TH1F*      fPDGCode ;		       //
    TH1F*      fGenCode ;		       // generator code
    TH1F*      fMom[2]  ;
    TH1F*      fTime    ;
    TH2F*      fYVsX    ;                // different VD's have different orientation
    TH2F*      fYVsZ    ;                // fill both hist's
    TH1F*      fPt      ;                // transverse mom
    TH1F*      fPp      ;                // momentum component parallel to the solenoid axis
    TH1F*      fTanTh   ;		       // tan (pitch angle)
    TH1F*      fEKin    ;
    TH2F*      fCosThVsMom ;	       // cos (pitch angle) vs Mom
    TH2F*      fCosThVsMomPV;		// for antiprotons
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
    EventHist_t*  fEvent      [kNEventHistSets      ];
    SimpHist_t*   fSimp       [kNSimpHistSets       ];
    VDetHist_t*   fVDet       [kNVDetHistSets       ];
    SpmcHist_t*   fStepPointMC[kNStepPointMCHistSets];
  };
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
public:
					// pointers to the data blocks used
  TSimpBlock*           fSimpBlock;  
  TStepPointMCBlock*    fSpmcBlock;
  TStepPointMCBlock*    fVDetBlock;
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
  StnVDetData_t         fVDet[200];
  int                   fStageID;
  int                   fNSimp;
  int                   fStage;

  int                   fNSpmc;
  int                   fNSpmc150;

  SimpData_t            fSimData[kMaxNSimp];

  TStntuple*            fStnt;
  double                fWeight; // event weight, determined by the production cross section

					// antiproton-specific : in the production vertex
  double                fPbarCosThPV;
  double                fPbarMomPV;

//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TStnSpmcAnaModule(const char* name="StnSpmcAna", const char* title="StnSpmcAna");
  ~TStnSpmcAnaModule();
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
  int     BeginJob();
  int     BeginRun();
  int     Event   (int ientry);
  int     EndJob  ();
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void    BookEventHistograms        (EventHist_t* Hist, const char* Folder);
  void    BookSimpHistograms         (SimpHist_t*  Hist, const char* Folder);
  void    BookStepPointMCHistograms  (SpmcHist_t*  Hist, const char* Folder);
  void    BookVDetHistograms         (VDetHist_t*  Hist, const char* Folder);

  void    FillEventHistograms        (EventHist_t* Hist);
  void    FillSimpHistograms         (SimpHist_t*  Hist, TSimParticle* Simp, SimpData_t* SimpData, double Weight = 1.);
  void    FillStepPointMCHistograms  (SpmcHist_t*  Hist, TStepPointMC* Step, SpmcData_t* SpmcData, double Weight = 1.);
  void    FillVDetHistograms         (VDetHist_t*  Hist, TStepPointMC* Step,                       double Weight = 1.);

  void    BookHistograms();
  void    FillHistograms();


  void    Debug();
//-----------------------------------------------------------------------------
// test
//-----------------------------------------------------------------------------
  void    Test001();

  ClassDef(TStnSpmcAnaModule,0)
};

#endif
