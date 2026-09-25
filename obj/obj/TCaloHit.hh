//-----------------------------------------------------------------------------
//  2026-09-09 P.Murat: Mu2e calorimeter hit data
//-----------------------------------------------------------------------------
#ifndef TCaloHit_hh
#define TCaloHit_hh

#include "TObject.h"
#include "TBuffer.h"

namespace mu2e {
  class CaloHit;
}

class TCaloHit : public TObject {
public:
  int            fCid;                  // crystal ID
  int            fNSipms;               // (number of R/O channels used, 1 or 2) || (n)digis) << 8
  int            fCrdIndex[2];          // == added in V2 == index of the CaloRecoDigi in the original reco list
  float          fTime;                 // 
  float          fEDep;                 //
  float          fSigT;                 // uncertainty on T.
  float          fSigE;                 // uncertainty on E
  
  mu2e::CaloHit* fOfflineCaloHit;       //! transient
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCaloHit();
  TCaloHit(int ID);
  virtual ~TCaloHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Cid       () const { return fCid;       }
  int     NSipms    () const { return (fNSipms     ) & 0xff; }
  int     NRecoDigis() const { return (fNSipms >> 8) & 0xff; }
  float   Time      () const { return fTime;      }
  float   EDep      () const { return fEDep;      }
  float   SigT      () const { return fSigT;      }
  float   SigE      () const { return fSigE;      }
  
  int     Disk      () const { return (fCid / 674); }
      
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "")       override;
  virtual void Print(Option_t* opt = "") const override;

  ClassDefOverride(TCaloHit,2)
};

#endif
