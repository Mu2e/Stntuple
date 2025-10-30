///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCosmicTrack_hh
#define TEvdCosmicTrack_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
// #include "TArc.h"
// #include "TEllipse.h"

namespace mu2e {
  class CosmicTrackSeed;
}

class TLine;

namespace stntuple {

class TEvdStrawHit;
class TEvdTrkStrawHit;

class TEvdCosmicTrack: public TObject {
public:
  
protected:
  int                           fNumber;
  const mu2e::CosmicTrackSeed*  fCTSeed;

  TObjArray*                    fListOfHits; // not owned, just cached pointers

  TLine*                        fLineXY;
  TLine*                        fLineZY;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdCosmicTrack();

  TEvdCosmicTrack(int Number, const mu2e::CosmicTrackSeed* CTSeed);

  virtual ~TEvdCosmicTrack();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int NHits() { return fListOfHits->GetEntriesFast(); }

  TEvdTrkStrawHit* Hit(int I) { 
    return (TEvdTrkStrawHit*) fListOfHits->UncheckedAt(I); 
  }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void AddHit(TObject* Hit) { fListOfHits->Add(Hit); }
//-----------------------------------------------------------------------------
// drawing functions
//-----------------------------------------------------------------------------
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive   (Int_t px, Int_t py) override;
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveVST(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveVRZ(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveCal(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint  (Option_t* Opt = "") override;
  virtual void  Clear  (Option_t* Opt = "") override;
  
  virtual void  PrintMe() const ;                           // **MENU**
  virtual void  Print  (Option_t* Opt = "") const override; // **MENU**

  ClassDef(stntuple::TEvdCosmicTrack,0)
};

}
#endif
