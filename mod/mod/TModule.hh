///////////////////////////////////////////////////////////////////////////////
//
//
#ifndef __Stntuple_mod_TModule_hh__
#define __Stntuple_mod_TModule_hh__

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"

#include <string>
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "Stntuple/obj/AbsEvent.hh"
// #include "Stntuple/print/TAnaDump.hh"

class TAnaRint;
class TAnaDump;

class TModule : public art::EDAnalyzer, public TNamed {

  enum { kNDebugBits = 100 };

public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>                       interactiveMode{Name("interactiveMode"), Comment("1: interactive mode"  ) };
      // fhicl::Sequence<std::string>           rootMacro      {Name("rootMacro"      ), Comment("good hit mask"        ) };
      fhicl::Atom<std::string>               rootMacro      {Name("rootMacro"      ), Comment("good hit mask"        ) };
      fhicl::Table<fhicl::ParameterSet>      debugBits      {Name("debugBits"      ), Comment("debug bits"           ) };
      //      fhicl::Table<TAnaDump::Config>         tanaDump       {Name("TAnaDump"       ), Comment("TAnaDumpParameters"   ) };
    };
					// there are some initializations which need 
					// to be done just once
  static int          fgInitialized;

  TFile*              fFile;

  int                 fDebugBit[kNDebugBits];		// flags for different debug options

  fhicl::ParameterSet fFclDebugBits;
  int                 fInteractiveMode;

  int                 fEventNumber;
  int                 fRunNumber;
  int                 fPassed;		// for filtering

  TAnaRint*           fAnaRint;
  TAnaDump*           fDump;

  const art::Run*     fRun;
					// provides for a possibility for any ROOT 
					// module to call (whenever necessary) a
					// loaded in interpreted function which 
					// name has to be defined at run time
					// it is a user responsibility to load a
					// macro with this function. Its signature:
					// 
					//       int function_name(int flag)
					// 
					// fFunction->Data() would return the 
					// `function_name'
  TString*            fFunction;
  std::string         _rootMacro;
					// each TModule adds a folder to gROOT
  TFolder*            fFolder;
//-----------------------------------------------------------------------------
// methods 
//-----------------------------------------------------------------------------
  TModule(const fhicl::ParameterSet&    pset, const char* Name);
  //   explicit TModule(const art::EDAnalyzer::Table<TModule::Config>& config, const char* Name);

  virtual      ~TModule();

  virtual int  beforeBeginJob();
  virtual void beginJob      ();
  virtual int  afterBeginJob ();

  virtual int  beforeBeginRun(const art::Run & _Run);
  virtual void beginRun      (const art::Run & _Run);
  virtual int  afterBeginRun (const art::Run & _Run);

  virtual int  beforeEvent   (const art::Event& Evt);
  virtual void analyze       (const art::Event& Evt);
  virtual int  afterEvent    (const art::Event& Evt);

  virtual int  beforeEndRun(const art::Run &  Rn);
  virtual void endRun      (const art::Run &  Rn);
  virtual int  afterEndRun (const art::Run &  Rn);

  virtual int  beforeEndJob();
  virtual void endJob      ();
  virtual int  afterEndJob ();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TFolder*     GetFolder() { return fFolder; }

  int          DebugBit(int I) { return fDebugBit[I]; }

  const char*  Function() { return (fFunction) ? fFunction->Data() : nullptr ; }

  int          ExecuteMacro(int Mode);
//-----------------------------------------------------------------------------
// the following helper methods allow to save 1 line per request, which in 
// case of 100's histograms booked is a non-negligible number
//-----------------------------------------------------------------------------
  void         AddHistogram(TObject* hist, const char* FolderName = "Hist");
  
  void         HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       const char* FolderName = "Hist");
  
  void         HBook2F(TH2F*& Hist, const char* Name, const char* Title,
		       Int_t Nx, Double_t XMin, Double_t XMax,
		       Int_t Ny, Double_t YMin, Double_t YMax,
		       const char* FolderName = "Hist");
  
  void         HProf (TProfile*& Hist, const char* Name, const char* Title,
		      Int_t Nx, Double_t XMin, Double_t XMax,
		      Double_t YMin, Double_t YMax,
		      const char* FolderName = "Hist");
  
  int          SaveFolder(TFolder* Folder, TDirectory* Dir);
  
};

#endif
