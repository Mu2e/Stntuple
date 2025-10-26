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

// #include "TrackerGeom/inc/Layer.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
// #include "TrackerGeom/inc/Sector.hh"

#include "Stntuple/gui/TEvdPanel.hh"
#include "Stntuple/gui/TEvdPlane.hh"
#include "Stntuple/gui/TEvdStraw.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TStnGeoManager.hh"

#include "Offline/TrackerGeom/inc/Panel.hh"
//#include "TTrackerGeom/inc/ZLayer.hh"
#include "Offline/DAQ/inc/TrkPanelMap_t.hh"

ClassImp(stntuple::TEvdPanel)

int stntuple::TEvdPanel::fgDebugLevel(0);

namespace stntuple {
//_____________________________________________________________________________
TEvdPanel::TEvdPanel(): TObject() {
  fVisible      = 0;
  fMnID         = -1;
}

//_____________________________________________________________________________
  TEvdPanel::TEvdPanel(int ID, const mu2e::Panel* Panel, TEvdPlane* Plane, const mu2e::Tracker* Tracker): TObject() {

  fID      = ID;
  fNLayers = Panel->nLayers();
  fPanel   = Panel;
  fVisible = 1;

  TStnGeoManager* gm = TStnGeoManager::Instance();
  const mu2e::TrkPanelMap::Row* tpm = gm->PanelMap(Plane->ID(),ID); // plane ID : 0-35 (???)
  if (tpm) fMnID              = tpm->mnid();
					// assume that the number of straws is the same
  int id;

  int nst = Panel->nStraws();
  fListOfStraws = new TObjArray(nst);

  mu2e::StrawId pid = Panel->id();
  if (TEvdPanel::fgDebugLevel != 0) {
    printf("--- station, plane, panel: %3i %3i %3i\n",
           pid.getStation(),pid.getPlane(),pid.getPanel());
  } 
  for (uint16_t ist=0; ist<nst; ist++) {
    const mu2e::Straw* straw = &Panel->getStraw(ist);

    id           = straw->id().asUint16();
    int ill      = straw->id().getLayer();
    int iss      = straw->id().getStraw();

    if (TEvdPanel::fgDebugLevel & 0x1) {
      printf(" layer, straw, z : %3i %3i %10.3f\n",ill,iss,straw->getMidPoint().z());
    }

    TEvdStraw* evd_straw = new TEvdStraw(id,straw,this,Tracker);
    fListOfStraws->Add(evd_straw);
  }

  if (TEvdPanel::fgDebugLevel & 0x2) {
//-----------------------------------------------------------------------------
// print rotation matrix and offset ?
//-----------------------------------------------------------------------------
    printf("origin: %10.4f  %10.4f  %10.4f \n",
           Panel->origin().x(),Panel->origin().y(),Panel->origin().z());
      
    printf("udir: %10.4f  %10.4f  %10.4f \n",
           Panel->uDirection().x(),Panel->uDirection().y(),Panel->uDirection().z());
      
    printf("vdir: %10.4f  %10.4f  %10.4f \n",
           Panel->vDirection().x(),Panel->vDirection().y(),Panel->vDirection().z());
      
    printf("wdir: %10.4f  %10.4f  %10.4f \n",
           Panel->wDirection().x(),Panel->wDirection().y(),Panel->wDirection().z());
  }
//-----------------------------------------------------------------------------
// build the transformation matrix
//-----------------------------------------------------------------------------
  double phiy   = atan2(Panel->vDirection().y(),Panel->vDirection().x())*180./M_PI;
  double phix   = phiy-90;
  double phiz   = 0;
  double thetax = 90;
  double thetay = 90;
  double thetaz =  0;
  if (Panel->wDirection().z() < 0) {
    phix   = phix+180;
  }

  // if ((pid.getPlane() % 2) == 1) {
  //   if ((pid.getPanel() % 2) == 0) {
  //     phix   = phix+180;
  //     phiy   = phiy+180;
  //   }
  // }
  // else {
  //   if ((pid.getPanel() % 2) == 1) {
  //     phix   = phix+180;
  //     phiy   = phiy+180;
  //   }
  // }

  TGeoRotation* r = new TGeoRotation("a",thetax,phix,thetay,phiy,thetaz,phiz);

  fCombiTrans     = new TGeoCombiTrans(Form("trk_%02i_%i",pid.getPlane(),pid.getPanel()),0,0,0,r);
  
  TVector3* p0  = Straw( 0)->Pos();
  TVector3* p95 = Straw(95)->Pos();

  double xc = (p0->X()+p95->X())/2;
  double yc = (p0->Y()+p95->Y())/2;
  double zc = (p0->Z()+p95->Z())/2;

  fPos.SetXYZ(xc,yc,zc);
}

//_____________________________________________________________________________
TEvdPanel::~TEvdPanel() {
}

//-----------------------------------------------------------------------------
void TEvdPanel::Paint(Option_t* Option) {


  int view = TVisManager::Instance()->GetCurrentView()->Type();

  if      (view == TStnVisManager::kXY) PaintXY (Option);
  else if (view == TStnVisManager::kRZ) PaintRZ (Option);
  else {
    // what is the default?
    //    Warning("Paint",Form("Unknown option %s",option));
  }

  gPad->Modified();
}


//-----------------------------------------------------------------------------
void TEvdPanel::PaintXY(Option_t* Option) {
}

//-----------------------------------------------------------------------------
void TEvdPanel::PaintRZ(Option_t* option) {
  // draw straws
}

//-----------------------------------------------------------------------------
  void TEvdPanel::PaintVST(Option_t* Option) {
  // draw straws
  int nstraws = NStraws();
  for (int is=0; is<nstraws; is++) {
    TEvdStraw* straw  = Straw(is);
    straw->PaintVST(Option);
  }
}

//-----------------------------------------------------------------------------
  void TEvdPanel::PaintVRZ(Option_t* Option) {
  // draw straws
  int nstraws = NStraws();
  for (int is=0; is<nstraws; is++) {
    TEvdStraw* straw  = Straw(is);
    straw->PaintVRZ(Option);
  }
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

  //  Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveVST(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdPanel::DistancetoPrimitiveVRZ(Int_t px, Int_t py) {
  return 9999;
}

}
