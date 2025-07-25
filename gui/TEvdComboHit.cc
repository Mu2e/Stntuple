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

#include "Offline/RecoDataProducts/inc/ComboHit.hh"

#include "Stntuple/print/TAnaDump.hh"

#include "Stntuple/gui/TEvdComboHit.hh"
#include "Stntuple/gui/TStnVisManager.hh"


ClassImp(stntuple::TEvdComboHit)

namespace stntuple {
//-----------------------------------------------------------------------------
TEvdComboHit::TEvdComboHit() {
}

//-----------------------------------------------------------------------------
TEvdComboHit::TEvdComboHit(const mu2e::ComboHit*      Hit,
			   const mu2e::SimParticle*   Sim,
			   const mu2e::StrawGasStep*  Step,
			   int                        MotherPdgID,
			   float                      P,
			   float                      Pz):
  TObject(),
  fHit(Hit),
  fPos(Hit->pos().x(),Hit->pos().y(),Hit->pos().z()),
  //  fDir(Hit->wdir().x(),Hit->wdir().y())
  fDir(Hit->uDir().x(),Hit->uDir().y())
{
  fSim         = Sim;
  fStep        = Step;
  if (Sim) {
    fSimID       = Sim->id().asInt();
    fPdgID       = Sim->pdgId();
  }
  else {
    fSimID       = -1;
    fPdgID       =  0;
  }
  fMotherPdgID = MotherPdgID;
  fP           = P;
  fPz          = Pz;

   // printf("TEvdComboHit::TEvdComboHit: StrawDidiMC::driftDistance disabled. ask Dave Brown\n");

  fSigW  = Hit->wireRes();
  fSigR  = 2.5 ;                // Hit->transRes();
//-----------------------------------------------------------------------------
// style and color
//-----------------------------------------------------------------------------
  int   style(0), color(0);
  float size(0.);

  if      (fPdgID == 11) {
    if    (fP      > 20  ) { style = 20; size = 0.8; color = kRed; }
    else                   { style = 20; size = 0.4; color = kRed+2; }
  }
  else if (fPdgID ==  -11) { style = 20; size = 0.8; color = kBlue;   }
  else if (fPdgID ==   13) { style = 20; size = 0.8; color = kGreen+2;}
  else if (fPdgID ==  -13) { style = 20; size = 0.8; color = kGreen-2;}
  else if (fPdgID ==  211) { style = 20; size = 0.8; color = kOrange-3;}
  else if (fPdgID == -211) { style = 20; size = 0.8; color = kOrange+7; }
  else if (fPdgID == 2212) { style = 20; size = 1.0; color = kBlue+2; }
  else                     { style = 20; size = 1.0; color = kBlack;}

  fColor = color;
//-----------------------------------------------------------------------------
// define lines for XY view
//-----------------------------------------------------------------------------
  fLineW.SetX1(fPos.X()-fDir.X()*fSigW);
  fLineW.SetY1(fPos.Y()-fDir.Y()*fSigW);
  fLineW.SetX2(fPos.X()+fDir.X()*fSigW);
  fLineW.SetY2(fPos.Y()+fDir.Y()*fSigW);
  fLineW.SetLineColor(color);

  fLineR.SetX1(fPos.X()+fDir.Y()*fSigR);
  fLineR.SetY1(fPos.Y()-fDir.X()*fSigR);
  fLineR.SetX2(fPos.X()-fDir.Y()*fSigR);
  fLineR.SetY2(fPos.Y()+fDir.X()*fSigR);
  fLineR.SetLineColor(color);

  //  const CLHEP::Hep3Vector* mid_point = &fStraw->GetStraw()->getMidPoint();

  // double rdrift(2.);
//-----------------------------------------------------------------------------
// RZ view - for straw hits
//-----------------------------------------------------------------------------
  // fEllipse.SetX1(mid_point->z());
  // fEllipse.SetY1(mid_point->perp());
  // fEllipse.SetR1(rdrift);
  // fEllipse.SetR2(rdrift);
  // fEllipse.SetFillStyle(3003);
  // fEllipse.SetFillColor(kBlue+2);
  // fEllipse.SetLineColor(kBlue+2);

  // fXYMarker.SetX(fX);
  // fXYMarker.SetY(fY);
  // fXYMarker.SetMarkerStyle(style);
  // fXYMarker.SetMarkerSize (size );
  // fXYMarker.SetMarkerColor(color);
//-----------------------------------------------------------------------------
// TZ view
//-----------------------------------------------------------------------------
  fTZMarker.SetX(fPos.Z());
  fTZMarker.SetY(Hit->correctedTime());
  fTZMarker.SetMarkerStyle(style);
  fTZMarker.SetMarkerSize (size );
  fTZMarker.SetMarkerColor(color);
//-----------------------------------------------------------------------------
// PhiZ view
//-----------------------------------------------------------------------------
  fPhiZMarker.SetX(fPos.Z());
  fPhiZMarker.SetY(fPos.Phi());
  fPhiZMarker.SetMarkerStyle(style);
  fPhiZMarker.SetMarkerSize (size );
  fPhiZMarker.SetMarkerColor(color);
}

//-----------------------------------------------------------------------------
TEvdComboHit::~TEvdComboHit() {
}

// //-----------------------------------------------------------------------------
// void TEvdComboHit::Paint(Option_t* Option) {
//   const char oname[] = "TEvdComboHit::Paint";

//   //  int   iv;

//   int view = TVisManager::Instance()->GetCurrentView()->Type();

//   if      (view == TStnView::kXY ) PaintXY (Option);
//   else if (view == TStnView::kTZ ) PaintTZ (Option);
//   else if (view == TStnView::kCal) PaintCal(Option);
//   else {
//     printf("[%s] >>> ERROR: unknown view: %i, DO NOTHING\n",oname,view);
//   }

//   gPad->Modified();
// }

//_____________________________________________________________________________
void TEvdComboHit::PaintXY(Option_t* Option) {
  fLineW.Paint(Option);
  fLineR.Paint(Option);
}

//_____________________________________________________________________________
void TEvdComboHit::PaintTZ(Option_t* Option) {
  fTZMarker.Paint(Option);
}

//_____________________________________________________________________________
void TEvdComboHit::PaintPhiZ(Option_t* Option) {
  fPhiZMarker.Paint(Option);
}

//_____________________________________________________________________________
void TEvdComboHit::Print(Option_t* Option) const {

  TAnaDump* ad = TAnaDump::Instance();

  int flags = *((int*) &fHit->flag());

  TString opt = Option;
  opt.ToLower();
  if (opt == "") ad->printComboHit(fHit, fStep, "banner+data", -1, flags);
  else           ad->printComboHit(fHit, fStep, Option       , -1, flags);

  // printf("TEvdComboHit::Print simID=%6i pdg=%10i motherPdgID=%10i time=%10.3f Z=%10.3f edep=%10.6f P=%10.3f Pz=%10.3f\n",
  // 	 fSimID,fPdgID,fMotherPdgID,fHit->time(),fPos.Z(),fHit->energyDep(),fP,fPz);
}

// //-----------------------------------------------------------------------------
// void TEvdComboHit::PaintCal(Option_t* option) {
//   // nothing to draw...
// }


// //_____________________________________________________________________________
// Int_t TEvdComboHit::DistancetoPrimitive(Int_t px, Int_t py) {
//   return 9999;
// }

// //-----------------------------------------------------------------------------
// Int_t TEvdComboHit::DistancetoPrimitiveXY(Int_t px, Int_t py) {

//   Int_t dist = 9999;

//   static TVector3 global;
// //   static TVector3 local;

// //   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

//   global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

//   return dist;
// }

// //_____________________________________________________________________________
// Int_t TEvdComboHit::DistancetoPrimitiveTZ(Int_t px, Int_t py) {
//   return 9999;
// }

}
