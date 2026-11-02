#ifndef __daqana_obj_TStrTrack_hh__
#define __daqana_obj_TStrTrack_hh__

#include "TObject.h"
#include "TVector3.h"
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
  std::vector<mu2e::TrkStrawHitSeed>* fHitSeeds;  // written if fFormat & 0x2 ==NOT OWNED==
  std::vector<mu2e::TrkStraw>*        fStraws;    // written if fFormat & 0x4 ==NOT OWNED==
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::KalSeed*                fOfflineKs; //! defined only during ntuple making
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  TStrTrack();
  TStrTrack(int ID, int Format = 1);
  virtual ~TStrTrack();

  int     NDof () const { return fNDof;  }
  int     NHits() const { return fNHits; }

  float   Nx   () const { return fNx; }
  float   Ny   () const { return fNy; }
  float   Nz   () const { return fNz; }

  float   T0   () const { return fT0; }

  float   X0   () const { return fX0; }
  float   Y0   () const { return fY0; }
  float   Z0   () const { return fZ0; }

  void    GetPosAtT(float T, float* Pos) const ;

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

  ClassDefOverride(TStrTrack,2);
};

#endif
