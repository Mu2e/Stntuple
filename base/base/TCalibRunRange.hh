//-----------------------------------------------------------------------------
// GetName() returns strings of format 'subdetector'/'table' 
// (for example, "ces/algn")
//-----------------------------------------------------------------------------
#ifndef TCalibRunRange_hh
#define TCalibRunRange_hh

#include "TObject.h"
#include "TString.h"

class TCalibRunRange: public TObject {
public:
  int      fMinRun;
  int      fMaxRun;
  TString  fFilename;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TCalibRunRange();
  TCalibRunRange(int MinRun, int MaxRun, const char* File);
  virtual ~TCalibRunRange();

  int         GetMinRun  () { return fMinRun; }
  int         GetMaxRun  () { return fMaxRun; }
  const char* GetFilename() { return fFilename.Data(); }

  ClassDefOverride(TCalibRunRange,1)
};

#endif

