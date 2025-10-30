///////////////////////////////////////////////////////////////////////////////
//
//

#include <iostream>
#include <boost/algorithm/string.hpp>
#include <string>

#include "Stntuple/mod/TModule.hh"
#include "Stntuple/mod/TAnaRint.hh"
#include "Stntuple/print/TAnaDump.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "TString.h"
#include "TROOT.h"
#include "TInterpreter.h"

using namespace std;

TModule::TModule(): art::EDAnalyzer(fhicl::ParameterSet{}), TNamed("","") {
}

  //-----------------------------------------------------------------------------
TModule::TModule(const fhicl::ParameterSet& PSet, const fhicl::ParameterSet& TModulePSet, const char* Name):
  art::EDAnalyzer(PSet),
  TNamed(Name,Name)
{

  int    n;
  static char* dummy[100];

  fFile   = 0;
  fFolder = new TFolder(GetName(),GetName());

  memset(fDebugBit,0,kNDebugBits*sizeof(int));

  fFclDebugBits    = TModulePSet.get<fhicl::ParameterSet>("debugBits"      );
  fInteractiveMode = TModulePSet.get<int>                ("interactiveMode");
  _rootMacro       = TModulePSet.get<string>             ("rootMacro"      );

  fAnaRint         = TAnaRint::Instance(0,dummy);
//-----------------------------------------------------------------------------
// pass parameters to TAnaDump
//-----------------------------------------------------------------------------
  fhicl::ParameterSet tadPset = TModulePSet.get<fhicl::ParameterSet>("TAnaDump");
  fDump                       = TAnaDump::Instance(&tadPset);

  const char* key;
                                        // a flag is an integer!
  n = fFclDebugBits.get_names().size();
  for (int i=0; i<n; i++) {
    int index(0);
    key               = fFclDebugBits.get_names().at(i).data();
    sscanf(key,"bit%i",&index );
    fDebugBit[index]  = fFclDebugBits.get<int>(key);

    printf("... TModule: bit=%3i is set to %i \n",index,fDebugBit[index]);
  }

//-----------------------------------------------------------------------------
// plugin initialization - figure out the function name
//-----------------------------------------------------------------------------
  fFunction = nullptr;
  if (_rootMacro.length() > 0) {
    TInterpreter::EErrorCode rc;
    TInterpreter* cint = gROOT->GetInterpreter();

    cint->LoadMacro(_rootMacro.data(),&rc);
    if (rc != 0) {
      printf("TModule:beginJob ERROR : failed to load %s\n",_rootMacro.data());
    }
    else {
      vector<string> res, r2;
      
      boost::algorithm::split(res,_rootMacro, boost::is_any_of("/"));

      int n = res.size();
					// last is the macro itself, strip the path
      string ww = res[n-1];

      boost::algorithm::split(r2,ww, boost::is_any_of("."));

      fFunction = new TString(r2[0].data());
    }
  }
}

//-----------------------------------------------------------------------------
TModule::TModule(const fhicl::Table<TModule::Config>& config, const char* Name):
  art::EDAnalyzer(config.get_PSet()),
  TNamed(Name,Name),
  fFclDebugBits   (config().debugBits()      ),
  fInteractiveMode(config().interactiveMode()),
  _rootMacro      (config().rootMacro()      )
{

  int    n, index;
  static char* dummy[100];

  fFile   = 0;
  fFolder = new TFolder(GetName(),GetName());

  memset(fDebugBit,0,kNDebugBits*sizeof(int));

  fAnaRint         = TAnaRint::Instance(0,dummy);
//-----------------------------------------------------------------------------
// pass parameters to TAnaDump
//-----------------------------------------------------------------------------
  fhicl::ParameterSet tadPset = config().TAnaDump.get<fhicl::ParameterSet>();

  fDump = TAnaDump::Instance(&tadPset);

  const char* key;
                                        // a flag is an integer!
  n = fFclDebugBits.get_names().size();
  for (int i=0; i<n; i++) {
    key                = fFclDebugBits.get_names().at(i).data();
    sscanf(key,"bit%i" ,&index );

    fDebugBit[index]  = fFclDebugBits.get<int>(key);

    printf("... TModule: bit=%3i is set to %i \n",index,fDebugBit[index]);
  }

//-----------------------------------------------------------------------------
// plugin initialization - figure out the function name
//-----------------------------------------------------------------------------
  fFunction = nullptr;
  if (_rootMacro.length() > 0) {
    TInterpreter::EErrorCode rc;
    TInterpreter* cint = gROOT->GetInterpreter();

    cint->LoadMacro(_rootMacro.data(),&rc);
    if (rc != 0) {
      printf("MuHitDisplay:beginJob ERROR : failed to load %s\n",_rootMacro.data());
    }
    else {
      vector<string> res, r2;
   
      boost::algorithm::split(res,_rootMacro, boost::is_any_of("/"));

      int n = res.size();
					// last is the macro itself, strip the path
      string ww = res[n-1];

      boost::algorithm::split(r2,ww, boost::is_any_of("."));

      fFunction = new TString(r2[0].data());
    }
  }
};

//-----------------------------------------------------------------------------
TModule::~TModule() {
  if (fFile) {
    delete fFile;
    fFile  = NULL;
  }
}


//______________________________________________________________________________
int TModule::beforeBeginJob() {
  return 0;
}

//-----------------------------------------------------------------------------
void TModule::beginJob() {
}


//______________________________________________________________________________
int TModule::afterBeginJob() {
  return 0;
}


//-----------------------------------------------------------------------------
int TModule::beforeBeginRun(const art::Run& aRun) {
  return 0;
}
    

//-----------------------------------------------------------------------------
void TModule::beginRun(const art::Run &  Rn) {
  ExecuteMacro(0);
}

//-----------------------------------------------------------------------------
int TModule::afterBeginRun(const art::Run& aRun) {
  return 0;
}
    
//______________________________________________________________________________
int TModule::beforeEvent(const AbsEvent& event) {
  return 0;
}

//-----------------------------------------------------------------------------
void TModule::analyze  (const art::Event& Evt) {
//-----------------------------------------------------------------------------
// if a ROOT macro is defined, execute it
//-----------------------------------------------------------------------------
  ExecuteMacro(1);

  if (fInteractiveMode == 1) {
    fDump->SetEvent(&Evt);
    fAnaRint->SetInteractiveMode(fInteractiveMode);
    fAnaRint->Rint()->Run(true);
                                        // provision for switching the interactive mode OFF

    fAnaRint->GetInteractiveMode(fInteractiveMode);
  }
}

//______________________________________________________________________________
int TModule::afterEvent(const AbsEvent& event) {
  return 0;
}

//______________________________________________________________________________
int TModule::beforeEndRun(const art::Run& _Run) {
  return 0;
}

//-----------------------------------------------------------------------------
// if fInteractiveMode=2, stop at the interactive prompt only in the end of run
// however the macro is executed for each event
//-----------------------------------------------------------------------------
void TModule::endRun(const art::Run &  Rn) {
  
  ExecuteMacro(2);

  if (fInteractiveMode == 2) {
					// at this point the event is not defined
    fDump->SetEvent(nullptr);
    fAnaRint->SetInteractiveMode(fInteractiveMode);
    fAnaRint->Rint()->Run(true);
                                        // provision for switching the interactive mode OFF

    fAnaRint->GetInteractiveMode(fInteractiveMode);
  }
}


//______________________________________________________________________________
int TModule::afterEndRun(const art::Run& _Run) {
  return 0;
}

//______________________________________________________________________________
int TModule::beforeEndJob() {
  return 0;
}

//-----------------------------------------------------------------------------
void TModule::endJob  () {
}

//______________________________________________________________________________
int TModule::afterEndJob() {
  return 0;
}

//_____________________________________________________________________________
void     TModule::AddHistogram(TObject* hist, const char* FolderName) {
  TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
  fol->Add(hist);
}


//-----------------------------------------------------------------------------
// Function(int Mode, TModule* Mod)
//-----------------------------------------------------------------------------
int TModule::ExecuteMacro(int Mode) {
  char  par[200];

  if (fFunction == nullptr) return 0;

  TInterpreter* cint = gROOT->GetInterpreter();
  // sprintf(par,"%s(%i,(TModule*) %p)",Function(),Mode,this);
  // cint->ProcessLine(par);

  sprintf(par,"%i,(TModule*)%p",Mode,static_cast<void*>(this));
  cint->Execute(Function(),par);

  return 0;
}


//_____________________________________________________________________________
  void TModule::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
                           Int_t Nx, Double_t XMin, Double_t XMax,
                           const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TH1F(Name,Title,Nx,XMin,XMax);
    AddHistogram(Hist,FolderName);
  }

//_____________________________________________________________________________
  void TModule::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
                           Int_t Nx, Double_t XMin, Double_t XMax,
                           Int_t Ny, Double_t YMin, Double_t YMax,
                           const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
    AddHistogram(Hist,FolderName);
  }

//_____________________________________________________________________________
  void TModule::HProf(TProfile*& Hist, const char* Name, const char* Title,
                         Int_t Nx, Double_t XMin, Double_t XMax,
                         Double_t YMin, Double_t YMax,
                         const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
    AddHistogram(Hist,FolderName);
  }


  //_____________________________________________________________________________
  int  TModule::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

    //    TFolder*     fol;
    TDirectory*  dir;
    TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
    Dir->cd();
    dir = Dir->mkdir(Folder->GetName(),Folder->GetName());
    dir->cd();

    TIter  it(Folder->GetListOfFolders());
    while ((o = it.Next())) {
      if (strcmp(o->ClassName(),"TFolder") == 0) {
        SaveFolder((TFolder*) o, dir);
      }
      else if (! o->InheritsFrom("TStnModule")) {
        o->Write();
        //      gDirectory->GetListOfKeys()->Print();
      }
    }

    Dir->cd();
    return 0;
  }
