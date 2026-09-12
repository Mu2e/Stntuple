//-----------------------------------------------------------------------------
//  2014-01-26 P.Murat: Mu2e calorimeter hit data
//-----------------------------------------------------------------------------
#ifndef TCalHitData_hh
#define TCalHitData_hh

#include "TObject.h"
#include "TBuffer.h"

namespace mu2e {
  class CalHit;
}

class TCalHitData : public TObject {
public:
  int            fCid;                  // crystal ID
  int            fNSipms;               // number of R/O channels used, 1 or 2
  float          fTime;                 // 
  float          fEDep;                 //
  float          fSigT;                 // uncertainty on T
  float          fSigE;                 // uncertainty on E
  
  mu2e::CalHit* fOfflineCalHit;       //! transient
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCalHitData();
  virtual ~TCalHitData();
					// ****** initialization
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     Cid      () { return fCid;       }
  int     ID       () { return fCid;       }   // obsolete
  int     NChannels() { return fNSipms;    }   // obsolete
  int     NSipms   () { return fNSipms;    }
  float   Time     () { return fTime;      }
  float   Energy   () { return fEDep;      }   // obsolete
  float   EDep     () { return fEDep;      }
  float   SigT     () { return fSigT;      }
  float   SigE     () { return fSigE;      }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void Set(int ID, int NSipms, float Time,  float EDep) {
    fCid = ID; fNSipms = NSipms; fTime = Time; fEDep = EDep;
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "")       override;
  virtual void Print(Option_t* opt = "") const override;

  ClassDefOverride(TCalHitData,2)
};

#endif
