#ifndef TDuplicateFilterModule_hh
#define TDuplicateFilterModule_hh

#include "TObjArray.h"
#include "Stntuple/loop/TStnModule.hh"

class TStnRunRecord;

class TDuplicateFilterModule: public TStnModule {
  // everything is public, it is communism
public:

  TObjArray*      fListOfRunRecords;
  TStnRunRecord*  fCurrentRunRecord;
  int             fNDuplicateEvents;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDuplicateFilterModule(const char* name  = "DuplicateFilter", 
			 const char* title = "DuplicateFilter");
  ~TDuplicateFilterModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// setters ... SetOutputDir leaks memory. I know.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  virtual int       BeginJob       () override;
  virtual int       BeginRun       () override;
  virtual int       Event          (int ientry) override;
  virtual int       EndJob         () override;

  ClassDefOverride(TDuplicateFilterModule,0)
};
#endif
