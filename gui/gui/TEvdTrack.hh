///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdTrack_hh
#define TEvdTrack_hh

// #include "Gtypes.h"
#include "TObjArray.h"
// #include "TH1.h"
// #include "TPad.h"
// #include "TArc.h"
#include "TEllipse.h"
#include "TLine.h"

namespace mu2e {
  class KalSeed;
}

namespace stntuple {

class TEvdStrawHit;
class TEvdTrkStrawHit;

class TEvdTrack: public TObject {
public:
  
protected:
  int                   fNumber;
  const mu2e::KalSeed*  fKSeed;

  TObjArray*            fListOfHits;

  TEllipse*             fEllipse;
  TLine*                fLineXY;
  TLine*                fLineZY;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdTrack();

  TEvdTrack(int Number, const mu2e::KalSeed* KSeed);

  virtual ~TEvdTrack();
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
  virtual void  PaintVRZ (Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive   (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveVRZ(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveCal(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");
  virtual void  Print(Option_t* Opt = "") const ; // **MENU**

  ClassDef(stntuple::TEvdTrack,0)
};

}
#endif
