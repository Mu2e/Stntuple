#ifndef __daqana_obj_TStrTrack_hh__
#define __daqana_obj_TStrTrack_hh__

#include "TObject.h"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkStraw.hh"

namespace mu2e {
  class KalSeed;
}

class TStrTrack : public TObject {
public:
  int     fNHits;
  int     fNDof ;
  int     fTcIndex;     // index of the corresponding time cluster
  int     fFormat;      // low byte: format
  
  float   fT0   ;
  float   fChi2 ;
  float   fX0   ;
  float   fY0   ;
  float   fZ0   ;
  float   fNx   ;
  float   fNy   ;
  float   fNz   ;
//-----------------------------------------------------------------------------
// added in V2
//-----------------------------------------------------------------------------
  std::vector<mu2e::TrkStrawHitSeed>*  fTrkStrawHitSeed; // written if fFormat & 0x2
  std::vector<mu2e::TrkStraw>*         fTrkStraw;        // written if fFormat & 0x4

const mu2e::KalSeed*  fOfflineKs; //!
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TStrTrack();
  TStrTrack(int ID, int Format = 1);
  virtual ~TStrTrack();

  int     NDof () { return fNDof;  }
  int     NHits() { return fNHits; }

  const mu2e::KalSeed* OfflineKs() { return fOfflineKs; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(const char* Opt = "")       override;
  virtual void Print(const char* Opt = "") const override;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);

  ClassDefOverride(TStrTrack,1);
};

#endif
