//-----------------------------------------------------------------------------
// GetName() returns strings of format 'subdetector'/'table' 
// (for example, "ces/algn")
//-----------------------------------------------------------------------------
#ifndef TCalibTable_hh
#define TCalibTable_hh

#include "TNamed.h"
#include "TObjArray.h"
#include "TClonesArray.h"

class TCalibRunRange;

class TCalibTable : public TNamed {
public:
  TObjArray*  fListOfRunRanges;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
public:
  TCalibTable();
  TCalibTable(const char* Name, const char* Detector);

  virtual ~TCalibTable();

  int   NRunRanges() const { return fListOfRunRanges->GetEntriesFast(); }

  TObjArray*  GetListOfRunRanges() { return fListOfRunRanges; }

  TCalibRunRange* GetRunRange(int RunNumber);

  const char* GetFilename(int RunNumber);

  ClassDefOverride (TCalibTable,1)
};

#endif

