///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdHelixSeed_hh
#define TEvdHelixSeed_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TEllipse.h"

#include "art/Framework/Principal/Event.h"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"

#include "Stntuple/gui/TStnVisNode.hh"

namespace stntuple {

class TEvdStrawHit;
class TEvdTrkStrawHit;

class TEvdHelixSeed: public TObject {
public:
  
protected:
  int                     fNumber;
  const mu2e::HelixSeed*  fHelixSeed;
  TStnVisNode*            fVisNode;	// backward link to the note
  TObjArray*              fListOfHits;
  TEllipse*               fEllipse;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdHelixSeed();

  TEvdHelixSeed(int Number, const mu2e::HelixSeed* HSeed, TStnVisNode* VisNode);

  virtual ~TEvdHelixSeed();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int NHits() { return fListOfHits->GetEntriesFast(); }

  TEvdTrkStrawHit* Hit(int I) { 
    return (TEvdTrkStrawHit*) fListOfHits->UncheckedAt(I); 
  }

  double T0() { return fHelixSeed->t0().t0(); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void AddHit(TObject* Hit) { fListOfHits->Add(Hit); }
//-----------------------------------------------------------------------------
// drawing functions
//-----------------------------------------------------------------------------
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive   (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveCal(Int_t px, Int_t py);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");
  virtual void  Print(Option_t* Opt = "") const ; // **MENU**

  ClassDefOverride(stntuple::TEvdHelixSeed,0)
};

}
#endif
