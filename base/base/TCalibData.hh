//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
#ifndef TCalibData_hh
#define TCalibData_hh

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"

class TCalibData : public TObject {
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
public:
  TCalibData();
  virtual ~TCalibData();

  virtual int Init(TObject* Data) = 0;

  ClassDefOverride (TCalibData,1)
};

#endif

