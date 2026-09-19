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
  int     Cid       () { return fCid;       }
  int     NSipms    () { return (fNSipms     ) & 0xff; }
  int     NRecoDigis() { return (fNSipms >> 8) & 0xff; }
  float   Time      () { return fTime;      }
  float   EDep      () { return fEDep;      }
  float   SigT      () { return fSigT;      }
  float   SigE      () { return fSigE;      }
  
  int     Disk      () { return (fCid / 674); }
      
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
