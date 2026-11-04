//--------------------------------------------------------------------------
// Description:
//	Class THistModule: base class for ROOT histogramming modules
//
// Environment: Mu2e
//
// Author List: P.Murat
//
// Copyright Information: 
//   Copyright (C) 1999		CDF/Fermilab
//------------------------------------------------------------------------

#include <assert.h>
#include <climits>
#include <cstdio>
#include <regex>
#include <string>

#include "TInterpreter.h"
#include "TH1.h"
#include "TFolder.h"
#include "TSystem.h"
#include "TError.h"

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/mod/THistModule.hh"

// ClassImp(THistModule)
					// initialize static class members

int        THistModule::fgBufferSize   = 64000;  // in bytes
int        THistModule::fgMaxFileSize  =   350;  // in MBytes
TString    THistModule::fgFileName     = "";
int        THistModule::fgSplitLevel   = 99;
int        THistModule::fgCompressionLevel = 1;

int        THistModule::fgFileNumber   = 0;
int        THistModule::fgOpenNextFile = 0;
TObjArray* THistModule::fgModuleList   = 0;
TFile*     THistModule::fgFile         = 0;
TTree*     THistModule::fgTree         = 0;
int        THistModule::fgMakeSubdirs  = 1;
int        THistModule::fgRunMin       = INT_MAX;
int        THistModule::fgSubrunMin    = INT_MAX;

//-----------------------------------------------------------------------------
THistModule::THistModule() : TModule() {
}

//______________________________________________________________________________
THistModule::THistModule(fhicl::ParameterSet const& PSet, fhicl::ParameterSet const& THistModulePSet, const char* Name): 
  TModule(PSet,THistModulePSet.get<fhicl::ParameterSet>("TModule"),Name)
{
  fOldDir        = NULL;
  fHistogramList = new TObjArray(10);

  fgBufferSize       = THistModulePSet.get<int>        ("bufferSize"      ,fgBufferSize);
  fgMaxFileSize      = THistModulePSet.get<int>        ("maxFileSize"     ,fgMaxFileSize);
  fgSplitLevel       = THistModulePSet.get<int>        ("splitLevel"      ,fgSplitLevel);
  fgCompressionLevel = THistModulePSet.get<int>        ("compressionLevel",fgCompressionLevel);

  TString s = THistModulePSet.get<std::string>("histFileName"    ,fgFileName.Data()).data();
  // redefine only if non-default
  if (s != "") fgFileName = s;

  std::cout << "-- module name:" << Name << " fgFileName:" << fgFileName.Data() << std::endl;

  fHistogramList->SetName("HistogramList");
  fFolder->Add(fHistogramList);

  fFolder->AddFolder("Hist","Hist");

  if (fgModuleList == 0) {
    fgModuleList = new TObjArray(10);
//-----------------------------------------------------------------------------
// do not append the histograms to the directory list of objects to allow
// using short histogram names
//-----------------------------------------------------------------------------
    TH1::AddDirectory(kFALSE);
  }
  fgModuleList->Add(this);
}

//-----------------------------------------------------------------------------
THistModule::THistModule(const fhicl::Table<THistModule::Config>& config, const char* Name) : 
  TModule                  (config.get_PSet(), Name)
  //  TModule(config.get_PSet(),Name)
{
  fOldDir        = NULL;
  fHistogramList = new TObjArray(10);

  fgBufferSize          = config().bufferSize     ();
  fgMaxFileSize         = config().maxFileSize     ();
  fgFileName            = config().histFileName    ();
  fgSplitLevel          = config().splitLevel      ();
  fgCompressionLevel    = config().compressionLevel();

  fHistogramList->SetName("HistogramList");
  fFolder->Add(fHistogramList);

  fFolder->AddFolder("Hist","Hist");

  if (fgModuleList == 0) {
    fgModuleList = new TObjArray(10);
    //    framework()->actions()->append(new THistModuleAction());
//-----------------------------------------------------------------------------
// do not append the histograms to the directory list of objects to allow
// using short histogram names
//-----------------------------------------------------------------------------
    TH1::AddDirectory(kFALSE);
  }
  fgModuleList->Add(this);
}

//______________________________________________________________________________
THistModule::~THistModule() {

//-----------------------------------------------------------------------------
//  close the file in the destructor such that all THistModules would
//  have their beforeEndJob entry points executed
//-----------------------------------------------------------------------------

  if (fgFile) {
    CloseAndStubFile(fgFile);
    fgFile = 0;
  }
  if (fHistogramList) {
    delete fHistogramList;
    fHistogramList = 0;
  }

  if (fgModuleList) {
    delete fgModuleList;
    fgModuleList = 0;
  }
}


//______________________________________________________________________________
void THistModule::ResetFileRange() {
  fgRunMin    = INT_MAX;
  fgSubrunMin = INT_MAX;
}

//______________________________________________________________________________
// lowest (run, subrun) pair of the events written into the current file
void THistModule::AddToFileRange(int Run, int Subrun) {
  if ((Run < fgRunMin) or ((Run == fgRunMin) and (Subrun < fgSubrunMin))) {
    fgRunMin    = Run;
    fgSubrunMin = Subrun;
  }
}

//______________________________________________________________________________
// close the file and rename it: the <run>_<subrun> field of the configured name
// (e.g. nts.mu2e.<dsid>.<cfg>.<run>_<subrun>.stn) is replaced with the lowest run
// and subrun of the events written into this file; the file counter suffix is kept
void THistModule::CloseAndStubFile(TFile* File) {
  if (File == nullptr) return;

  TString openName = File->GetName();
  File->Write();
  delete File;                                   // the TFile destructor also closes the file

  if (fgRunMin == INT_MAX) return;               // no events written into this file

  std::string base = fgFileName.Data();
  std::regex  stub("\\.[0-9]{6}_[0-9]{6}\\.");
  if (not std::regex_search(base, stub)) {
    ::Error("CloseAndStubFile","no <run>_<subrun> field in %s, file name not changed",fgFileName.Data());
    return;
  }

  char repl[32];
  snprintf(repl,sizeof(repl),".%06d_%06d.",fgRunMin,fgSubrunMin);
  TString newName = std::regex_replace(base, stub, repl, std::regex_constants::format_first_only).c_str();

  TString suffix = openName;                     // file counter suffix, e.g. "_2", empty for the first file
  suffix.Remove(0,fgFileName.Length());
  newName += suffix;

  if (newName != openName) {
    if (gSystem->Rename(openName.Data(),newName.Data()) != 0) {
      ::Error("CloseAndStubFile","failed to rename %s to %s",openName.Data(),newName.Data());
    }
  }
}

//______________________________________________________________________________
int THistModule::OpenNewFile(const char* Filename) {
  // open next file
  int rc = 0;
  fgFile = new TFile(Filename,"RECREATE");
  if (! fgFile) {
    Error("beginJob","an attempxt to open a new ROOT file %s failed",
	  fgFileName.Data());
    rc = -1;
  }
  return rc;
}


//______________________________________________________________________________
int THistModule::beforeBeginJob() {

					// return code
  int rc  =  0;
				// give more time to define TModule::fName
  fDirName = GetName();

  if ((fgFileName != "") && (fgFile == 0)) {

					// create a new ROOT File 
    rc = OpenNewFile(fgFileName.Data());
  }

  if (fgFile) {
					// file opened, don't forget to make
					// new directory before booking 
    if (fgMakeSubdirs) {
      fOldDir = gDirectory;
      fgFile->mkdir(fDirName);
      fgFile->cd(fDirName);
    }
  }

  rc = TModule::beforeBeginJob();

  return rc;
}


//______________________________________________________________________________
int THistModule::afterBeginJob() {

  if (fgMakeSubdirs && fgFile) {
    fOldDir->cd();
  }

  int rc = TModule::afterBeginJob();
  return rc;
}


//______________________________________________________________________________
int THistModule::beforeBeginRun(const art::Run& _Run) {
  return TModule::beforeBeginRun(_Run);
}
    

//______________________________________________________________________________
int THistModule::afterBeginRun(const art::Run& _Run) {
  return TModule::afterBeginRun(_Run);
}


//______________________________________________________________________________
int THistModule::beforeEvent(const AbsEvent& event) {
  if (fgFile) {
					// don't forget to change the directory
					// before starting filling histograms
    fOldDir = gDirectory;

    if (fgOpenNextFile) {
//-----------------------------------------------------------------------------
// current file will be closed by FillStntupleModule, this event will be 
// written into the next one, so write out the folder, 
// then reset the histograms
//-----------------------------------------------------------------------------
      fgFile->cd();
      fFolder->Write();
      fOldDir->cd();

      TObjArray* list = (TObjArray*) fFolder->FindObject("Hist");
      if (list) {
	int nhist = list->GetEntriesFast();
	for (int i=0; i<nhist; i++) {
	  TH1* h = (TH1*) list->UncheckedAt(i);
	  h->Reset();
	}
      }
    }
//-----------------------------------------------------------------------------
//  I suspect fgMakeSubdirs is always 1 now
//-----------------------------------------------------------------------------
    if (fgMakeSubdirs) {
      fgFile->cd(fDirName);
    }
  }

  int rc = TModule::beforeEvent(event);
  return rc;
}
    

//______________________________________________________________________________
int THistModule::afterEvent(const AbsEvent& _Event) {
  // need to cd() back 

  TModule::afterEvent(_Event);

  if (fgFile) {
    if (fgMakeSubdirs) {
      fOldDir->cd();
    }
  }
  int rc = TModule::afterEvent(_Event);
  return rc;
}
    

//______________________________________________________________________________
int THistModule::beforeEndRun(const art::Run& _Run) {
  int rc = TModule::beforeEndRun(_Run);
  return rc;
}
    

//______________________________________________________________________________
int THistModule::afterEndRun(const art::Run& _Run) {
  int rc = TModule::afterEndRun(_Run);
  return rc;
}
    

//______________________________________________________________________________
int THistModule::beforeEndJob() {
  if (fgFile) {
					// don't forget to change the directory
					// before starting filling histograms
    fOldDir = gDirectory;

    fgFile->cd();
    fFolder->Write();
    fOldDir->cd();
  }

  int rc = TModule::beforeEndJob();
  return rc;
}
    

//-----------------------------------------------------------------------------
int THistModule::afterEndJob() {
					// return code
  return TModule::afterEndJob();
}
    

//_____________________________________________________________________________
int THistModule::SetFileName(const char* nm) {
  if (fgFile) {
    printf(">>> ERROR: ROOT histogram file is already open. Renaming is disabled.\n");
  }
  else {
    if (nm && (fgFileName  != "")) {
      printf(">>> WARNING : %s is setting a new name for ROOT histogram file : %s\n",
	     GetName(), nm);
    }
    fgFileName = nm;
  }
  return 0;
}

