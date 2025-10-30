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

#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TStnGeoManager.hh"

#include "Offline/TrackerGeom/inc/Plane.hh"

// #include "Offline/DAQ/inc/TrkPanelMap_t.hh"


ClassImp(stntuple::TEvdPlane)

namespace stntuple {
//_____________________________________________________________________________
TEvdPlane::TEvdPlane(): TObject() {
  fID           = -1;
  fListOfPanels = NULL;
  fNPanels      = 0;
  fVisible      = 0;
  fProductionID = -1;
}

//_____________________________________________________________________________
  TEvdPlane::TEvdPlane(int ID, const mu2e::Plane* Plane, TEvdStation* Station, const mu2e::Tracker* Tracker): TObject() {

  TEvdPanel*  evd_panel;

  fID      = ID;
  fStation = Station;
  fPlane   = Plane;
  fNPanels  = Plane->nPanels();
  fVisible  = 1;

  fListOfPanels = new TObjArray(fNPanels);

  TStnGeoManager* gm          = TStnGeoManager::Instance();
  const mu2e::TrkPanelMap::Row* tpm = gm->PanelMap(ID,0);           // first panel of the plane
  if (tpm) fProductionID      = tpm->ppid();

  for (int i=0; i<fNPanels; i++) {
    const mu2e::Panel* panel = &fPlane->getPanel(i);

    evd_panel = new TEvdPanel(i,panel,this,Tracker);

    fListOfPanels->Add(evd_panel);
  }
}

//_____________________________________________________________________________
TEvdPlane::~TEvdPlane() {
  fListOfPanels->Delete();
  delete fListOfPanels;
}

//-----------------------------------------------------------------------------
void TEvdPlane::Paint(Option_t* option) {
  // paints one disk (.. or vane, in the past), i.e. section

				// parse option list

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
void TEvdPlane::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdPlane::PaintRZ(Option_t* option) {
}

//_____________________________________________________________________________
Int_t TEvdPlane::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdPlane::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdPlane::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

}
