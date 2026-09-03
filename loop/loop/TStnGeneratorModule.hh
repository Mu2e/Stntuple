#ifndef STNTUPLE_TStnGeneratorModule_hh
#define STNTUPLE_TStnGeneratorModule_hh

#include "Stntuple/loop/TStnInputModule.hh"

class TPrimaryInteraction;
class TG3Generator;
class TMcGenEvent;

class TStnGeneratorModule: public TStnInputModule {
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
protected:
  TMcGenEvent*          fMcEvent;
  TString               fBeam;
  TString               fTarget;
  TString               fFrame;
  Float_t               fEnergy;
  TObjArray*            fListOfGenerators;
  int                   fRunNumber;           // run number to be generated
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TStnGeneratorModule(TG3Generator* Gen = NULL); 
  virtual ~TStnGeneratorModule();

  TMcGenEvent* GetMcEvent() { return fMcEvent; }
//-----------------------------------------------------------------------------
//  accessors
//-----------------------------------------------------------------------------
  TG3Generator* Generator(const char* Name) { 
    return (TG3Generator*) fListOfGenerators->FindObject(Name); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetFrame (const char* Name  ) { fFrame  = Name;   }
  void SetBeam  (const char* Name  ) { fBeam   = Name;   }
  void SetTarget(const char* Name  ) { fTarget = Name;   }
  void SetEnergy(double      Energy) { fEnergy = Energy; }
  void SetRunNumber(int RunNumber) { fRunNumber = RunNumber; }

  void AddGenerator(TG3Generator* Gen) { 
    fListOfGenerators->Add((TObject*) Gen); 
  }
//-----------------------------------------------------------------------------
//  overloaded methods of TStnInputModule
//-----------------------------------------------------------------------------
  virtual int      NextEvent (Int_t I);
  virtual Double_t GetEntries() { return DBL_MAX; }
//-----------------------------------------------------------------------------
//  overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int BeginJob    ();
  virtual int BeginRun    ();
  virtual int Event       (Int_t I);
  virtual int EndRun      ();
  virtual int EndJob      ();

  virtual TStnNode* GetNode(const char* BranchName, const char* ClassName);

  virtual Int_t RegisterInputBranches(TStnEvent* Event);

  ClassDefOverride(TStnGeneratorModule,0) // MC event generator module
};

#endif
