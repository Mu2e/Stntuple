#ifndef __daqana_obj_TCaloRecoDigi_hh__
#define __daqana_obj_TCaloRecoDigi_hh__

#include "TObject.h"

namespace mu2e {
  class CaloRecoDigi;
}

class TCaloRecoDigi : public TObject {
public:
  int    fSipmID;
  int    fNdf;
  int    fPileup;
  float  fEDep;
  float  fSigE;
  float  fTime;
  float  fSigT;
  float  fChi2;
  const mu2e::CaloRecoDigi*  fOfflineCrd; //!
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TCaloRecoDigi();
  TCaloRecoDigi(int ID);
  virtual ~TCaloRecoDigi();

  int     Ndf   () { return fNdf;   }
  int     SipmID() { return fSipmID;}
  float   Time  () { return fTime;  }

  const mu2e::CaloRecoDigi* OfflineCrd() { return fOfflineCrd; }
  
  virtual void Clear(const char* Opt = "")       override;
  virtual void Print(const char* Opt = "") const override;

  ClassDefOverride(TCaloRecoDigi,1);
};

#endif
