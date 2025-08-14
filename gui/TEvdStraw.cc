///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
///////////////////////////////////////////////////////////////////////////////
#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLine.h"
#include "TArc.h"
#include "TArrow.h"
#include "TMath.h"
#include "TBox.h"
#include "TObjArray.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TEvdTimeCluster.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(stntuple::TEvdStraw)

namespace stntuple {
//_____________________________________________________________________________
TEvdStraw::TEvdStraw(): TObject() {
  fListOfHits = new TObjArray();
}

//-----------------------------------------------------------------------------
// drawing straws makes sense only on RZ view, in any other view it is just 
// a waste of time ... and screen space
//-----------------------------------------------------------------------------
  TEvdStraw::TEvdStraw(int Index, const mu2e::Straw* Straw, TEvdPanel* Panel,
                       const mu2e::Tracker* Tracker):  TObject()  {
  fIndex      = Index;
  fStraw      = Straw;
  fPanel      = Panel;
  fHalfLength = fStraw->halfLength();

  const CLHEP::Hep3Vector& pos = fStraw->getMidPoint  ();
  const CLHEP::Hep3Vector& dir = fStraw->wireDirection();

  fPos.SetXYZ(pos.x(),pos.y(),pos.z());
  fDir.SetXYZ(dir.x(),dir.y(),dir.z());

  double z     = fPos.Z();
  double rwire = pos.perp();                      // radial position of the wire
  double r     = Tracker->strawProperties()._strawOuterRadius;

  fArc = new TArc(z,rwire,r);

  double x1,y1, x2,y2;

  x1 = fPos.X()-fHalfLength*fDir.X();
  y1 = fPos.Y()-fHalfLength*fDir.Y();

  x2 = fPos.X()+fHalfLength*fDir.X();
  y2 = fPos.Y()+fHalfLength*fDir.Y();

  fWire = new TLine(x1,y1,x2,y2);

  fListOfHits = new TObjArray(5);
  fListOfHits->SetOwner(kFALSE);
}

//_____________________________________________________________________________
TEvdStraw::~TEvdStraw() {

  fListOfHits->Clear();
  delete fListOfHits;

  delete fArc;
}

//-----------------------------------------------------------------------------
void TEvdStraw::Paint(Option_t* option) {
  // paints one disk (.. or vane, in the past), i.e. section

  int view = TVisManager::Instance()->GetCurrentView()->Type();

  if      (view == TStnVisManager::kXY) PaintXY (option);
  else if (view == TStnVisManager::kRZ) PaintRZ (option);
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  gPad->Modified();
}


//_____________________________________________________________________________
void TEvdStraw::PaintXY(Option_t* Option) {
}


//_____________________________________________________________________________
void TEvdStraw::PaintRZ(Option_t* Option) {
  // draw straw

  int    nhits, color(0), style(0);
  double xstraw, ystraw, x1, x2, y1, y2;

  TStnVisManager* vm = TStnVisManager::Instance();
  float tmin = vm->TMin(); 
  float tmax = vm->TMax();
  stntuple::TEvdTimeCluster* etcl = vm->SelectedTimeCluster();
  if (etcl) {
    tmin = etcl->TMin(); // FIXME!
    tmax = etcl->TMax(); // FIXME!
  }
//-----------------------------------------------------------------------------
// do not even attempt to paint straws outside the visible range
//-----------------------------------------------------------------------------
  gPad->GetRange(x1,y1,x2,y2);

  xstraw = fArc->GetX1();
  ystraw = fArc->GetY1();

  if ((xstraw < x1) || (xstraw > x2) || (ystraw < y1) || (ystraw > y2)) return;

  nhits = fListOfHits->GetEntriesFast();
//-----------------------------------------------------------------------------
// the hit drift time is unknown, so just change the color of the straw circle
//-----------------------------------------------------------------------------
  fArc->SetLineColor(1);
  fArc->SetLineWidth(1);

  for (int i=0; i<nhits; ++i) {
    stntuple::TEvdStrawHit* evd_sh = Hit(i);
    const mu2e::ComboHit*   sch    = evd_sh->StrawHit();
    float t = sch->correctedTime();
    if ((t >= tmin) and (t <= tmax)) {    
      int ok = 1;
      if (etcl and vm->DisplayOnlyTCHits()) { 
        ok = etcl->TCHit(sch->index());
      }
      if (ok) {
        color = kRed;
        fArc->SetLineColor(color+1);
        fArc->SetLineWidth(2);
      }
    }
  }

  fArc->SetFillStyle(style);
  fArc->Paint(Option);

  if (vm->DisplayStrawDigiMC()) {
    for (int i=0; i<nhits; i++) {
      Hit(i)->PaintRZ(Option);
    }
  }
}

//_____________________________________________________________________________
void TEvdStraw::PaintVST(Option_t* Option) {
  double x1,y1, x2,y2;
  gPad->GetRange(x1,y1,x2,y2);

  fWire->Paint(Option);
}

//-----------------------------------------------------------------------------
// RZ view of the straw is the same, just keep the namning scheme
//-----------------------------------------------------------------------------
void TEvdStraw::PaintVRZ(Option_t* Option) {
  PaintRZ(Option);
}

//_____________________________________________________________________________
Int_t TEvdStraw::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdStraw::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdStraw::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdStraw::DistancetoPrimitiveVST(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdStraw::DistancetoPrimitiveVRZ(Int_t px, Int_t py) {
  return 9999;
}


//-----------------------------------------------------------------------------
void TEvdStraw::Clear(const char* Opt) {
  fListOfHits->Clear();
}

//-----------------------------------------------------------------------------
void TEvdStraw::Print(const char* Opt) const {
  printf("[TEvdStraw::Print] ID:%5i NHits: 2i\n",fNHits);
}

}
