#ifndef __daqana_obj_TStrTrack_hh__
#define __daqana_obj_TStrTrack_hh__

#include "TObject.h"

namespace mu2e {
  class KalSeed;
}

class TStrTrack : public TObject {
public:
  int     fNHits;
  int     fNDof ;
  int     fTcIndex;     // index of the corresponding time cluster
  
  float   fT0   ;
  float   fChi2 ;
  float   fX0   ;
  float   fY0   ;
  float   fZ0   ;
  float   fNx   ;
  float   fNy   ;
  float   fNz   ;
  const mu2e::KalSeed*    fOfflineKs; //!
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TStrTrack();
  TStrTrack(int ID);
  virtual ~TStrTrack();

  int     NDof () { return fNDof;  }
  int     NHits() { return fNHits; }

  const mu2e::KalSeed* OfflineKs() { return fOfflineKs; }
  
  virtual void Clear(const char* Opt = "")       override;
  virtual void Print(const char* Opt = "") const override;

  ClassDefOverride(TStrTrack,1);
};

#endif
