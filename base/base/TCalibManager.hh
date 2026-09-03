//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
#ifndef TCalibManager_hh
#define TCalibManager_hh

#include <cstdio>

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

class TCalibTable;
class TCalibRunRange;

class TCalibManager : public TObject {
public:
  static TCalibManager* fgInstance;
  TString               fDirectory;
  TString               fPass;
  TObjArray*            fListOfTables; // ! list of cached tables (internal)
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
public:
  TCalibManager();
  virtual ~TCalibManager();

  virtual  TCalibTable*    GetTable   (const char* Detector, const char* Table);

  virtual  TCalibRunRange* GetRunRange(const char* Detector, 
				       const char* Table,
				       int         RunNumber);

  static TCalibManager* Instance() {
    return (fgInstance) ? fgInstance : new TCalibManager();
  }

  int  SetPass(const char* Pass, Int_t Force = 0);

  ClassDefOverride (TCalibManager,1)
};

#endif

