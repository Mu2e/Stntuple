///////////////////////////////////////////////////////////////////////////////
// May 04 2013 P.Murat
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

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// #include "Offline/GeometryService/inc/GeometryService.hh"
// #include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdStation.hh"
#include "Stntuple/gui/TEvdTracker.hh"
#include "Stntuple/gui/TStnVisManager.hh"

ClassImp(stntuple::TEvdTracker)

namespace stntuple {
//_____________________________________________________________________________
  TEvdTracker::TEvdTracker(const mu2e::Tracker* Tracker): TNamed("tracker","Mu2e tracker") {

  fTracker = Tracker;

  if (Tracker == NULL) {
    printf(">>> TEvdTracker::TEvdTracker ERROR: Tracker = NULL\n");
  }

  fNStations      = mu2e::StrawId::_nstations; // fTracker->nStations();
  fListOfStations = new TObjArray(fNStations);

  for (int i=0; i<fNStations; i++) {
    TEvdStation* s = new TEvdStation(i,Tracker);
    fListOfStations->Add(s);
  }
}

//-----------------------------------------------------------------------------
// need a destructor...
//-----------------------------------------------------------------------------
TEvdTracker::~TEvdTracker() {
  if (fListOfStations) {
    fListOfStations->Delete();
    delete fListOfStations;
  }
}

//-----------------------------------------------------------------------------
void TEvdTracker::Paint(Option_t* option) {

  int view = TVisManager::Instance()->GetCurrentView()->Type();

  if      (view == TStnVisManager::kXY ) PaintXY (option);
  else if (view == TStnVisManager::kRZ ) PaintRZ (option);
  else if (view == TStnVisManager::kVST) PaintVST(option);
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  gPad->Modified();
}


//_____________________________________________________________________________
void TEvdTracker::PaintXY(Option_t* Option) {
}



//_____________________________________________________________________________
void TEvdTracker::PaintRZ(Option_t* option) {
  // draw tracker

  TEvdStation*   station;
  TEvdPlane*     plane;
  //  TEvdFace*      face;
  TEvdPanel*     panel;
  TEvdStraw*     straw;

  int            nplanes, npanels, nstraws;

  for (int ist=0; ist<fNStations; ist++) {
    station = Station(ist);
    nplanes = station->NPlanes();
    for (int ipl=0; ipl<nplanes; ipl++) {
      plane   = station->Plane(ipl);
      npanels = plane->NPanels();
      for (int ipanel=0; ipanel<npanels; ipanel++) {
	panel   = plane->Panel(ipanel);
        nstraws = panel->NStraws();
        for (int is=0; is<nstraws; is++) {
          straw  = panel->Straw(is);
          straw->PaintRZ();
	}
      }
    }
  }
}

//-----------------------------------------------------------------------------
// only Panel needs PaintVST
//-----------------------------------------------------------------------------
void TEvdTracker::PaintVST(Option_t* Option) {
  printf("TEvdTracker::%s is not implemented yet.\n",__func__);

  TEvdStation*   station;
  TEvdPlane*     plane;
  TEvdPanel*     panel;

  int            nplanes, npanels;

  for (int ist=0; ist<fNStations; ist++) {
    station = Station(ist);
    if (station->Visible() == 0) continue;
    
    nplanes = station->NPlanes();
    for (int ipl=0; ipl<nplanes; ipl++) {
      plane   = station->Plane(ipl);
      if (plane->Visible() == 0) continue;
      npanels = plane->NPanels();
      for (int ipanel=0; ipanel<npanels; ipanel++) {
	panel   = plane->Panel(ipanel);
        if (panel->Visible() == 0) continue;
	panel->PaintVST(Option);
      }
    }
  }
  
}

//_____________________________________________________________________________
Int_t TEvdTracker::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdTracker::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdTracker::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdTracker::DistancetoPrimitiveVST(Int_t px, Int_t py) {
  return 9999;
}

}
