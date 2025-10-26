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
#include "TEllipse.h"
#include "TObjArray.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "art/Framework/Principal/Handle.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Offline/TrackerGeom/inc/Straw.hh"
//#include "TrackerConditions/inc/Types.hh"

#include "Stntuple/print/TAnaDump.hh"

#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdStrawHit.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(stntuple::TEvdStrawHit)

namespace stntuple {
//-----------------------------------------------------------------------------
TEvdStrawHit::TEvdStrawHit() {
}

//-----------------------------------------------------------------------------
TEvdStrawHit::TEvdStrawHit(const mu2e::ComboHit*    Hit,
			   TEvdStraw*               Straw,
			   const mu2e::StrawDigiMC* Sdmc,
			   double X, double Y, double Z, 
			   double                   Wx,
			   double                   Wy,
			   double                   SigW,
			   double                   SigR,
			   int                      Mask, 
			   int                      Color): 
  TObject(),
  fHit(Hit),
  fSdmc(Sdmc),
  fStraw(Straw),
  fPos(X,Y,Z),
  fDir(Wx,Wy),
  fEllipse()
 {
   // printf("TEvdStrawHit::TEvdStrawHit: StrawDidiMC::driftDistance disabled. ask Dave Brown\n");

  fSigW  = SigW;
  fSigR  = SigR;
  fMask  = Mask;
//-----------------------------------------------------------------------------
// style and color
//-----------------------------------------------------------------------------
  fColor = Color;
//-----------------------------------------------------------------------------
// define marker and lines
//-----------------------------------------------------------------------------
  fMarker.SetX(fPos.X());
  fMarker.SetY(fPos.Y());
  fMarker.SetMarkerSize(0.4);
  fMarker.SetMarkerColor(Color);

  fLineW.SetX1(fPos.X()-fDir.X()*fSigW);
  fLineW.SetY1(fPos.Y()-fDir.Y()*fSigW);
  fLineW.SetX2(fPos.X()+fDir.X()*fSigW);
  fLineW.SetY2(fPos.Y()+fDir.Y()*fSigW);
  fLineW.SetLineColor(Color);

  fLineR.SetX1(fPos.X()+fDir.Y()*fSigR);
  fLineR.SetY1(fPos.Y()-fDir.X()*fSigR);
  fLineR.SetX2(fPos.X()-fDir.Y()*fSigR);
  fLineR.SetY2(fPos.Y()+fDir.X()*fSigR);
  fLineR.SetLineColor(Color);

  const CLHEP::Hep3Vector* mid_point = &fStraw->GetStraw()->getMidPoint();
  
  double rdrift(0.);
  if (fSdmc) {
    if (fSdmc->wireEndTime(mu2e::StrawEnd::cal) < fSdmc->wireEndTime(mu2e::StrawEnd::hv)) {
      rdrift = 5.; // fSdmc->driftDistance(mu2e::StrawEnd::cal);
    }
    else {
      rdrift = 5.; // fSdmc->driftDistance(mu2e::StrawEnd::hv);
    }
  }
      
  fEllipse.SetX1(mid_point->z());
  fEllipse.SetY1(mid_point->perp());
  fEllipse.SetR1(rdrift);
  fEllipse.SetR2(rdrift);
  fEllipse.SetFillStyle(3003);
  fEllipse.SetFillColor(kBlue-10);
  fEllipse.SetLineColor(kBlue-10);
}

//-----------------------------------------------------------------------------
TEvdStrawHit::~TEvdStrawHit() {
}

//-----------------------------------------------------------------------------
void TEvdStrawHit::Paint(Option_t* Option) {
  const char oname[] = "TEvdStrawHit::Paint";

  //  int   iv;

  int view = TVisManager::Instance()->GetCurrentView()->Type();

  if      (view == TStnVisManager::kXY ) PaintXY (Option);
  else if (view == TStnVisManager::kRZ ) PaintRZ (Option);
  else if (view == TStnVisManager::kCal) PaintCal(Option);
  else {
    printf("[%s] >>> ERROR: unknown view: %i, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//_____________________________________________________________________________
void TEvdStrawHit::PaintXY(Option_t* Option) {
  fMarker.Paint(Option);
  fLineW.Paint(Option);
  fLineR.Paint(Option);
}

//_____________________________________________________________________________
void TEvdStrawHit::PaintRZ(Option_t* Option) {
  if (fHit->flag().hasAllProperties(mu2e::StrawHitFlag::active)) {
    fEllipse.SetFillColor(kBlue-10);
  }
  else {
    fEllipse.SetFillColor(kRed-10);
  }
  fEllipse.SetFillStyle(3004);
  fEllipse.Paint(Option);
}

//_____________________________________________________________________________
void TEvdStrawHit::PaintCal(Option_t* option) {
  // nothing to draw...
}


// //_____________________________________________________________________________
// Int_t TEvdStrawHit::DistancetoPrimitive(Int_t px, Int_t py) {
//   return 9999;
// }

//_____________________________________________________________________________
Int_t TEvdStrawHit::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdStrawHit::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//-----------------------------------------------------------------------------
void TEvdStrawHit::Print(Option_t* Option) const {

  TAnaDump* ad = TAnaDump::Instance();

  const mu2e::StrawGasStep* step(nullptr);

  if (fSdmc) { 
    if (fSdmc->wireEndTime(mu2e::StrawEnd::cal) < fSdmc->wireEndTime(mu2e::StrawEnd::hv)) {
      step = fSdmc->strawGasStep(mu2e::StrawEnd::cal).get();
    }
    else {
      step = fSdmc->strawGasStep(mu2e::StrawEnd::hv ).get();
    }
  }

  int flags = *((int*) &fHit->flag());

  TString opt = Option;
  opt.ToLower();
  if (opt == "") ad->printComboHit(fHit, step, "banner+data", -1, flags);
  else           ad->printComboHit(fHit, step, Option       , -1, flags);
}

//-----------------------------------------------------------------------------
void TEvdStrawHit::PrintMe() const {
  Print("");
}

}
