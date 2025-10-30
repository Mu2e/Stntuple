//--------------------------------------------------------------------------
// File and Version Information: RootHistModule.hh,v 1.0
//
// Description:
//	Class THistModule: base class for ROOT histogramming modules
//
// Environment: CDF Run II
//
// Author List: P.Murat
//
// Copyright Information: 
//   Copyright (C) 1999		CDF/Fermilab
//------------------------------------------------------------------------

#ifndef THistModule_HH
#define THistModule_HH

#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/mod/TModule.hh"

class TTree;

class THistModule : public TModule {
//------------------------------------------------------------------------------
//  static data members
//------------------------------------------------------------------------------
public:
#ifndef __CLING__
    struct Config { 
      using Name    = fhicl::Name; 
      using Comment = fhicl::Comment;

      fhicl::Atom<int>                     bufferSize      {Name("bufferSize"      ),Comment("buffer size"       ) };
      fhicl::Atom<int>                     maxFileSize     {Name("maxFileSize"     ),Comment("max file size"     ) };
      // fhicl::Sequence<TString>       histFileName    {Name("histFileName"    ),Comment("hist file name"    ) };
      fhicl::Atom<TString>                 histFileName    {Name("histFileName"    ),Comment("hist file name"    ) };
      fhicl::Atom<int>                      splitLevel      {Name("splitLevel"      ),Comment("split level"       ) };
      fhicl::Atom<int>                      compressionLevel{Name("compressionLevel"),Comment("compression level" ) };
    };
#endif
protected:
					// there are some initializations 
					// which need to be done just once
  static TString    fgFileName;
  static TFile*     fgFile;
  static TObjArray* fgModuleList;
  static TTree*     fgTree;
  static int        fgMakeSubdirs;
  static int        fgMaxFileSize;
  static int        fgFileNumber;
  static int        fgOpenNextFile;
  static int        fgSplitLevel;
  static int        fgBufferSize;
  static int        fgCompressionLevel;
//------------------------------------------------------------------------------
//  data members of the module
//------------------------------------------------------------------------------
					// name of the directory in a ROOT file
					// associated with the module
  TString           fDirName;
					// list of histograms/ntuples owned by
					// the module
  TObjArray*        fHistogramList;
					// cache for cd() command
  TDirectory*       fOldDir;
//------------------------------------------------------------------------------
//  methods of the class
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  THistModule();
  
#ifndef __CLING__  
  explicit THistModule(const fhicl::ParameterSet&   PSet  ,
                       const fhicl::ParameterSet&   THistModulePSet,
                       const char*                  Name);
  explicit THistModule(const fhicl::Table<THistModule::Config>& Config, const char* Name);
#endif
  ~THistModule( );
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int        SplitLevel      () { return fgSplitLevel;   }
  int        CompressionLevel() { return fgCompressionLevel; }
  int        BufferSize      () { return fgBufferSize;       }
  TFile*     File            () { return fgFile;         }
  TObjArray* HistogramList   () { return fHistogramList; }


  static:: TObjArray* GetListOfModules() { return fgModuleList; }

//-----------------------------------------------------------------------------
// other methods
// define name of the output ROOT file, this can be done once per job
//-----------------------------------------------------------------------------
  int        SetFileName    (const char* Filename);
  int        OpenNewFile    (const char* Filename);

					// ****** overloaded methods of the 
					// base class

					// use beginJob to book histograms
					// use event to fill histograms

					// ****** methods called by the action
					// controller

  virtual int beforeBeginJob();
  virtual int afterBeginJob ();
  virtual int beforeBeginRun(const art::Run& aRun );
  virtual int afterBeginRun (const art::Run& aRun );
  virtual int beforeEvent   (const AbsEvent& event);
  virtual int afterEvent    (const AbsEvent& event);
  virtual int beforeEndRun  (const art::Run& aRun );
  virtual int afterEndRun   (const art::Run& aRun );
  virtual int beforeEndJob  ();
  virtual int afterEndJob   ();

  //  ClassDef(THistModule,0)
};

#endif
