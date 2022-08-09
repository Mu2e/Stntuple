///////////////////////////////////////////////////////////////////////////////
// May 04 2022 P.Murat
// 
// in 'XY' mode draw calorimeter clusters as circles with different colors 
// in 'Cal' mode draw every detail...
//
// BaBar interface:
// ----------------
//      r     = fabs(1./om);
//      phi0  = Trk->helix(0.).phi0();
//      x0    =  -(1/om+d0)*sin(phi0);
//      y0    =   (1/om+d0)*cos(phi0);
//
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
#include "TPolyLine.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"


#include "BTrk/KalmanTrack/KalRep.hh"

#include "art/Framework/Principal/Handle.h"

#include "BTrk/KalmanTrack/KalRep.hh"

#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/HelixTraj.hh"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"

#include "Offline/TrackerGeom/inc/Tracker.hh"

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

#include "Stntuple/gui/TEvdSimParticle.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TEvdTrkStrawHit.hh"
#include "Stntuple/base/TObjHandle.hh"

#include "CLHEP/Vector/ThreeVector.h"

ClassImp(stntuple::TEvdSimParticle)

namespace stntuple {

//-----------------------------------------------------------------------------
TEvdSimParticle::TEvdSimParticle(): TObject() {
  fListOfHits  = NULL;
  fEllipse     = new TEllipse();
  fParticlePDG = nullptr;
  fStep        = nullptr;
  fSimp        = nullptr;
}

//-----------------------------------------------------------------------------
// pointer to track is just cached
// Step - at the trackerffront
//-----------------------------------------------------------------------------
  TEvdSimParticle::TEvdSimParticle(int Number, const mu2e::SimParticle* Simp, const mu2e::StepPointMC* Step): TObject() {
  fNumber = Number;
  fSimp   = Simp;
  fStep   = Step;

  TDatabasePDG* pdb = TDatabasePDG::Instance();

  fParticlePDG = pdb->GetParticle(fSimp->pdgId());

  fListOfHits = new TObjArray();

  fEllipse = new TEllipse();

  double r, phi0, x0, y0, xc, yc, p, pt, q;

  x0 = (fStep->position().x()+3904.); // in mm
  y0 = fStep->position().y();

  p  = fStep->momentum().mag();
  pt = fStep->momentum().perp();
  r  = pt/2.9979*10;                          //    10^10/c, in mm

  phi0 =  fStep->momentum().phi();
  q    = fParticlePDG->Charge()/3;	// returned is the charge in units of 1/3 

    
  xc   =  x0+r*q*sin(phi0);
  yc   =  y0-r*q*cos(phi0);

  // // printf("[MuHitDispla::printHelixParams] d0 = %5.3f r = %5.3f phi0 = %5.3f x0 = %5.3f y0 = %5.3f\n",
  // // 	 d0, r, phi0, x0, y0);
   
  fEllipse->SetR1(r);
  fEllipse->SetR2(r);
  fEllipse->SetX1(xc);
  fEllipse->SetY1(yc);
  fEllipse->SetPhimin(0);
  fEllipse->SetPhimax(2*M_PI*180);
  fEllipse->SetTheta(0);

  int pdg_id = fSimp->pdgId();

  int color;

  if      (pdg_id == 11) {
    if    (p       > 20  ) { color = kRed;    }
    else                   { color = kRed+2;  }
  }
  else if (pdg_id ==  -11) { color = kBlue;   } 
  else if (pdg_id ==   13) { color = kGreen+2;} 
  else if (pdg_id ==  -13) { color = kGreen-2;} 
  else if (pdg_id == 2212) { color = kBlue+2; } 
  else                     { color = kMagenta;} 

  fEllipse->SetFillStyle(0);
  fEllipse->SetFillColor(0);
					// figure out the color
  fEllipse->SetLineColor(color);
  // fEllipse->SetFillStyle(3001);		// make it transparent
    
}

//-----------------------------------------------------------------------------
TEvdSimParticle::~TEvdSimParticle() {
  delete fEllipse;

  if (fListOfHits) {
    fListOfHits->Delete();
    delete fListOfHits;
  }
}

//-----------------------------------------------------------------------------
void TEvdSimParticle::Paint(Option_t* Option) {
  const char oname[] = "TEvdSimParticle::Paint";

  int view = TVisManager::Instance()->GetCurrentView()->Type();

  if      (view == TStnView::kXY ) PaintXY (Option);
  else if (view == TStnView::kRZ ) PaintRZ (Option);
  else if (view == TStnView::kCal) {
//-----------------------------------------------------------------------------
// calorimeter-specific view: do not draw tracks
//-----------------------------------------------------------------------------
  }
  else {
    printf("[%s] >>> ERROR: unknown view: %i, DO NOTHING\n",oname,view);
  }

  gPad->Modified();
}

//-----------------------------------------------------------------------------
// to display the reconstructed track in XY, use its parameters in the middle 
// of the tracker, at s=0
//-----------------------------------------------------------------------------
void TEvdSimParticle::PaintXY(Option_t* Option) {

  // fEllipse->PaintEllipse(x0,y0,r,r,0,2*M_PI*180,0);
  fEllipse->Paint();
}

//-----------------------------------------------------------------------------
void TEvdSimParticle::PaintRZ(Option_t* Option) {

  //  double            flen, zwire[2], ds, rdrift, zt[4], rt[4], zw, rw;
  double            zwire[2], rdrift, zt[4], rt[4], zw, rw;
  CLHEP::Hep3Vector tdir;
  HepPoint          tpos;
  TPolyLine         pline;
  int               nplanes, /*npanels,*/ nl;

  mu2e::GeomHandle<mu2e::Tracker> handle;

  const mu2e::Tracker* tracker = handle.get();

  //  const mu2e::Straw  *hstraw, *s, *straw[2];		// first straw
  const mu2e::Straw  *hstraw, *straw[2];		// first straw
  
  nplanes = tracker->nPlanes();
//-----------------------------------------------------------------------------
// first display track hits - active and not 
//-----------------------------------------------------------------------------
  const mu2e::TrkStrawHit  *hit;
  //   const TrkHitVector*       hits = &fKrep->hitVector();

  int nhits = 0; // NHits();
  for (int i=0; i<nhits; i++) {
    hit = nullptr; // Hit(i)->TrkStrawHit();
    rdrift = hit->driftRadius();

    hstraw = &hit->straw();
    zw     = hstraw->getMidPoint().z();
    rw     = hstraw->getMidPoint().perp();

    fEllipse->SetX1(zw);
    fEllipse->SetY1(rw);
    fEllipse->SetR1(0);
    fEllipse->SetR2(rdrift);
    fEllipse->SetFillStyle(3003);

    if (hit->isActive()) fEllipse->SetFillColor(kRed);
    else                 fEllipse->SetFillColor(kBlue+3);

    fEllipse->PaintEllipse(zw,rw,rdrift,0,0,2*M_PI*180,0);
  }
//-----------------------------------------------------------------------------
// now draw the trajectory in Z-R(local) space - device is a station
//-----------------------------------------------------------------------------
//  HelixTraj trkHel(fKrep->helix(0).params(),fKrep->helix(0).covariance());

  for (int iplane=0; iplane<nplanes; iplane++) {
    const mu2e::Plane* plane = &tracker->getPlane(iplane);
//-----------------------------------------------------------------------------
// a plane is made of 2 'faces' or 6 panels, panels 0 and 1 on different faces
//-----------------------------------------------------------------------------
//    npanels = plane->nPanels();
//-----------------------------------------------------------------------------
// 3 panels in the same plane have the same Z and do not overlap - 
// no need to duplicate; panels 0 and 1 are at different Z
//-----------------------------------------------------------------------------
    for (int iface=0; iface<2; iface++) {
      const mu2e::Panel* panel = &plane->getPanel(iface);
      nl       = panel->nLayers();
					// deal with the compiler warnings
      zwire[0] = -1.e6;
      zwire[1] = -1.e6;

      for (int il=0; il<nl; il++) {
//-----------------------------------------------------------------------------
// assume all wires in a layer have the same Z, extrapolate track to the layer
//-----------------------------------------------------------------------------
	straw[il] = &panel->getStraw(il);
	zwire[il] = straw[il]->getMidPoint().z();
      }
					// order locally in Z
      if (zwire[0] > zwire[1]) {
	zt[1] = zwire[1];
	zt[2] = zwire[0];
	//	flip  = 1.;
      }
      else {
	zt[1] = zwire[0];
	zt[2] = zwire[1];
	//	flip  = -1.;
      }
					// z-step between the layers, want two more points
      zt[0] = zt[1]-3.;
      zt[3] = zt[2]+3.;

      // const CLHEP::Hep3Vector* wd;
      // double r;
      
//       for (int ipoint=0; ipoint<4; ipoint++) {
// //-----------------------------------------------------------------------------
// // estimate flen using helix
// //-----------------------------------------------------------------------------
// 	flen   = trkHel.zFlight(zt[ipoint]);
// 	fKrep->traj().getInfo(flen,tpos,tdir);
// //-----------------------------------------------------------------------------
// // try to extrapolate helix to a given Z a bit more accurately 
// //-----------------------------------------------------------------------------
// 	ds     = (zt[ipoint]-tpos.z())/tdir.z();
// 	fKrep->traj().getInfo(flen+ds,tpos,tdir);
// //-----------------------------------------------------------------------------
// // for each face, loop over 3 panels in it
// //-----------------------------------------------------------------------------
// 	rt[ipoint] = -1.e6;
// 	for (int ipp=0; ipp<3; ipp++) {
// 	  const mu2e::Panel* pp = &plane->getPanel(2*ipp+iface);
// 	  s  = &pp->getStraw(0);
// 	  wd = &s->getDirection();
// 	  r  = (tpos.x()*wd->y()-tpos.y()*wd->x()); // *flip;
// 	  if (r > rt[ipoint]) rt[ipoint] = r;
// 	}
//       }
      pline.SetLineWidth(2);
      pline.PaintPolyLine(4,zt,rt);
    }
  }

}

//_____________________________________________________________________________
Int_t TEvdSimParticle::DistancetoPrimitive(Int_t px, Int_t py) {
  return 9999;
}

//_____________________________________________________________________________
Int_t TEvdSimParticle::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t dist = 9999;

  static TVector3 global;
//   static TVector3 local;

//   Double_t    dx1, dx2, dy1, dy2, dx_min, dy_min, dr;

  global.SetXYZ(gPad->AbsPixeltoX(px),gPad->AbsPixeltoY(py),0);

  return dist;
}

//_____________________________________________________________________________
Int_t TEvdSimParticle::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}


//-----------------------------------------------------------------------------
void TEvdSimParticle::Clear(Option_t* Opt) {
  //  SetLineColor(1);
}


//-----------------------------------------------------------------------------
void TEvdSimParticle::Print(Option_t* Opt) const {

//   TObjHandle*                  h;
//   const mu2e::CaloCrystalHit*  hit;

//   printf (" X0 = %10.3f Y0 = %10.3f  E = %10.3f  njits = %5i\n",
// 	  X0(),Y0(),fEnergy,fNHits);

//   printf("----------------------------------------------------------------\n");
//   printf("CrystalID      Time   Energy    EnergyTot  NRoids               \n");
//   printf("----------------------------------------------------------------\n");
//   for (int i=0; i<fNHits; i++) {

//     h   = (TObjHandle*) fListOfHits->At(i);
//     hit = (const mu2e::CaloCrystalHit*) h->Object();

//     printf("%7i  %10.3f %10.3f %10.3f %5i\n",
// 	   hit->id(),
// 	   hit->time(),
// 	   hit->energyDep(),
// 	   hit->energyDepTotal(),
// 	   hit->numberOfROIdsUsed());
//   }
}
}