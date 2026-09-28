#ifndef __daqana_obj_TCaloRecoDigi_hh__
#define __daqana_obj_TCaloRecoDigi_hh__

#include "TObject.h"
#include "TBuffer.h"

namespace mu2e {
  class CaloRecoDigi;
}

class TCaloRecoDigi : public TObject {
public:
  int    fSipmID;
  int    fNdf;
  int    fPileup;
  int    fCdIndex;      // == added in V2 == (calodigi index)
  float  fEDep;
  float  fSigE;
  float  fTime;
  float  fSigT;
  float  fChi2;
  float  fNext;                           //! first transient, marks end of record
  const mu2e::CaloRecoDigi*  fOfflineCrd; //!
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TCaloRecoDigi();
  TCaloRecoDigi(int ID);
  virtual ~TCaloRecoDigi();

  int     Ndf    () const { return fNdf;   }
  int     SipmID () const { return fSipmID;}
  int     CdIndex() const { return fCdIndex;}
  float   Time   () const { return fTime;  }

  const mu2e::CaloRecoDigi* OfflineCrd() { return fOfflineCrd; }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(const char* Opt = "")       override;
  virtual void Print(const char* Opt = "") const override;

  ClassDefOverride(TCaloRecoDigi,2);
};

#endif
