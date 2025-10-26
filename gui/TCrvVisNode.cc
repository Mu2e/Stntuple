#include <cmath>

#include "TVirtualX.h"
#include "TPad.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TBox.h"
#include "TObjArray.h"
#include "TColor.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"

#include "Stntuple/gui/TCrvVisNode.hh"
#include "Stntuple/gui/TStnVisManager.hh"
#include "Stntuple/gui/TEvdCrvBar.hh"

#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorModule.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorBar.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"


ClassImp(TCrvVisNode)

//-----------------------------------------------------------------------------
TCrvVisNode::TCrvVisNode(const char* Name, int SectionID) : TStnVisNode(Name) {
//-----------------------------------------------------------------------------
// This gradient palette will always be the same
// but the time-scale is governed by the vis manager slider.
// That way, a time window can be selected, and the colors will
// go from blue-red (low-high) in the time window
//-----------------------------------------------------------------------------
  Double_t stop [] = { 0.00, 0.20, 0.40, 0.60, 0.70, 1.00 };
  Double_t r    [] = { 0.00, 0.00, 0.00, 0.97, 0.97, 0.10 };
  Double_t g    [] = { 0.97, 0.30, 0.40, 0.97, 0.00, 0.00 };
  Double_t b    [] = { 0.97, 0.97, 0.00, 0.00, 0.00, 0.00 };
  Int_t    FI      = TColor::CreateGradientColorTable(6, stop, r, g, b, 1000);

  for (int i = 0; i < 1000; i++) colorPalette[i] = FI + i;

  fMinPulsePEs = 10; // Default number for minimum PEs

  mu2e::GeomHandle< mu2e::CosmicRayShield > CRS;
  int nmodules, nlayers, nbars;

  TEvdCrvBar* bar;
  
  fSectionID = SectionID;
  
  fListOfBars = new TObjArray();

  int ns = CRS->getCRSScintillatorShields().size(); // ????
  int sh09( 9);
  int sh13(13); // //
  int sh14(14); // //
 
  switch (fSectionID) {
  case 0:
    // Right
    for (int shield = 0; shield <=5; shield++) { // Loop over all the shields in right, but skip short section
      if (shield >= ns) continue;
      nmodules = CRS->getCRSScintillatorShield(shield).nModules();
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
	nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
	for (int layer = 0; layer < nlayers; layer++) {
	  nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
	  for (int ib = 0; ib < nbars; ib++) { 
	    bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(ib), fSectionID);
	    fListOfBars->Add(bar);
	  }
	}
      }
      if (!shield) { // Skip short shield (1)
	shield++;
      }
    }
    break;

  case 1:
    // Left
    for (int shield = 6; shield <= 8; shield++) { // Loop over all the shields in left
      if (shield >= ns) continue;
      nmodules = CRS->getCRSScintillatorShield(shield).nModules();
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
	nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
	for (int layer = 0; layer < nlayers; layer++) { // Loop over all the layers in the module
	  nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
	  for (int ib=0; ib< nbars; ib++) {
	    bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(ib), fSectionID);
	    fListOfBars->Add(bar);
	  }
	}
      }
    }
    break;
    
  case 2:
//-----------------------------------------------------------------------------
// Top DS
//-----------------------------------------------------------------------------
    for (int shield = 10; shield <= 12; shield++) { // Loop over all the shields in Top DS
      if (shield >= ns) continue;
      nmodules = CRS->getCRSScintillatorShield(shield).nModules();
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
	nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
	for (int layer = 0; layer < nlayers; layer++) { // Loop over all the layers in the module
	  nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
	  for (int ib=0; ib<nbars; ib++) {
	    bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(ib), fSectionID);
	    fListOfBars->Add(bar);
	  }
	}
      }
    }
    break;

  case 3:
//-----------------------------------------------------------------------------
// Downstream
//-----------------------------------------------------------------------------
    if (sh13 < ns) {
      nmodules = CRS->getCRSScintillatorShield(sh13).nModules(); // Get downstream shield
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
        const mu2e::CRSScintillatorModule& mod = CRS->getCRSScintillatorShield(sh13).getModule(module);
        nlayers = mod.nLayers();
        for (int layer = 0; layer < nlayers; layer++) { // Loop over all the layers in the module
          nbars = mod.getLayer(layer).nBars();
          for (int ib = 0; ib < nbars; ib++) {
            bar = new TEvdCrvBar(mod.getLayer(layer).getBar(ib), fSectionID);
            fListOfBars->Add(bar);
          }
        }
      }
    }
//-----------------------------------------------------------------------------
// Downstream lower sections - added 05/04/2015
//-----------------------------------------------------------------------------
    for (int shield = 18; shield <= 20; shield++) {
      if (shield >= ns) continue;
      nmodules = CRS->getCRSScintillatorShield(shield).nModules();
      for (int module = 0; module < nmodules; module++) {
	nlayers = CRS->getCRSScintillatorShield(shield).getModule(module).nLayers();
	for (int layer = 0; layer < nlayers; layer++) {
	  nbars = CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).nBars();
	  for (int ib = 0; ib < nbars; ib++) {
	    bar = new TEvdCrvBar(CRS->getCRSScintillatorShield(shield).getModule(module).getLayer(layer).getBar(ib), fSectionID);
	    fListOfBars->Add(bar);
	  }
	}
      }
    }
    break;
    
  case 4:
    // Upstream
    if (sh14 < ns) {
      nmodules = CRS->getCRSScintillatorShield(sh14).nModules(); // Get usptream shield
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
        const mu2e::CRSScintillatorModule& mod = CRS->getCRSScintillatorShield(sh14).getModule(module);
        nlayers = mod.nLayers();
        for (int layer = 0; layer < nlayers; layer++) { // Loop over all the layers in the module
          nbars = mod.getLayer(layer).nBars();
          for (int ib = 0; ib < nbars; ib++) { // Loop over all the bars in the layer
            bar = new TEvdCrvBar(mod.getLayer(layer).getBar(ib), fSectionID);
            fListOfBars->Add(bar);
          }
        }
      }
    }
    break;

  case 5: // Cryo Upstream
  case 6: // Cryo Downstream
  case 7: // Cryo Top
    // Cryo hole shields are not displayed for now
    break;

  case 8:
    // Top TS
    if (sh09 < ns) {
      nmodules = CRS->getCRSScintillatorShield(sh09).nModules(); // Get Top TS shield
      for (int module = 0; module < nmodules; module++) { // Loop over all the modules in the shield
        const mu2e::CRSScintillatorModule& mod = CRS->getCRSScintillatorShield(sh09).getModule(module);
        nlayers = mod.nLayers();
        for (int layer = 0; layer < nlayers; layer++) { // Loop over all the layers in the module
          nbars = mod.getLayer(layer).nBars();
          for (int ib = 0; ib < nbars; ib++) {
            bar = new TEvdCrvBar(mod.getLayer(layer).getBar(ib), fSectionID);
            fListOfBars->Add(bar);
          }
        }
      }
    }
    break;

  default:
    Warning("TCrvVisNode", Form("Unknown SectionID %i", SectionID));
    break;
  }
}

//-----------------------------------------------------------------------------
TCrvVisNode::~TCrvVisNode() {
}

//-----------------------------------------------------------------------------
int TCrvVisNode::InitEvent() {	
  ftimeLow = 400;
  ftimeHigh = 1695;

  TEvdCrvBar*              evd_bar;
  mu2e::GeomHandle<mu2e::CosmicRayShield> CRS;

  Clear();

  if (fCrvRecoPulsesCollection) { //If we have a valid pointer to a collection (non-empty), then proceed to fill
    //Loop over the RecoPulses in the collection
    for (mu2e::CrvRecoPulseCollection::const_iterator icrpc = (*fCrvRecoPulsesCollection)->begin(), 
	   ecrpc = (*fCrvRecoPulsesCollection)->end(); icrpc != ecrpc; ++icrpc) {
      const mu2e::CRSScintillatorBar &CRVCounterBar = CRS->getBar(icrpc->GetScintillatorBarIndex());

      //If the bar that we are looking at it is not a bar in this CRV Section, skip it
      int sec = getCRVSection(CRVCounterBar.id().getShieldNumber());
      if (sec != fSectionID) continue;

      //Set the bar pointer to the bar with map bar index
      int ind = icrpc->GetScintillatorBarIndex().asInt();
      evd_bar = EvdBarWithIndex(ind);

      //	  const mu2e::CrvRecoPulse &crvRecoPulses = icrpc->second; //The set of pulses for this bar

      int SiPM = icrpc->GetSiPMNumber();
//-----------------------------------------------------------------------------
// Add the pulse to the bar.
// If the pulse is above the threshold and within the initial time window, color the sipm
//-----------------------------------------------------------------------------
      evd_bar->AddPulse(*icrpc, SiPM); 
      if ((icrpc->GetPEs()       > fMinPulsePEs) && 
	  (icrpc->GetPulseTime() > ftimeLow    ) && 
	  (icrpc->GetPulseTime() < ftimeHigh)  )   { 
	
	evd_bar->SetFillStyle(1001, SiPM);
	evd_bar->SetFillColor(colorPalette[(int) ((icrpc->GetPulseTime() - ftimeLow) / (ftimeHigh - ftimeLow) * 999)], SiPM);
      }
	   
      //Loop over each SiPM for this bar
      // for (unsigned int SiPM = 0; SiPM < 4; SiPM++)
      //   {
      //     const std::vector<mu2e::CrvRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPM);
      
      //     //Loop over single pulses
      //     for (unsigned int i = 0; i < pulseVector.size(); i++)
      // 	{
      // 	  const mu2e::CrvRecoPulse &pulse = pulseVector[i];
      // 	  evd_bar->AddPulse(pulse, SiPM); // Add the pulse to the bar
      
      // 	  if ((pulse._PEs > fMinPulsePEs) && (pulse._pulseTime > ftimeLow) && (pulse._pulseTime < ftimeHigh)) 
      // If the pulse is above the threshold and within the initial time window, color the sipm
      // 	    {						
      // 	      evd_bar->SetFillStyle(1001, SiPM);
      // 	      evd_bar->SetFillColor(colorPalette[(int) ((pulse._pulseTime - ftimeLow) / (ftimeHigh - ftimeLow) * 999)], SiPM);
      // 	    }
      // 	} // Loop over single pulses
      //   } // Loop over SiPMs
    }
  }
  //  printf("Finished TCrvVisNode::InitEvent() for section %i \n", fSectionID);
  
  return 0;
}

//-----------------------------------------------------------------------------
void TCrvVisNode::UpdateEvent() {
  printf("Updating event... TCrvVisNode::UpdateEvent() for section %i \n", fSectionID);
  
  int nbars = NBars();
  for (int ibar = 0; ibar<nbars; ibar++) {
    TEvdCrvBar* bar = EvdBar(ibar);
    bar->SetThreshold(fMinPulsePEs);
    bar->SetTimeLow  (ftimeLow);
    bar->SetTimeHigh (ftimeHigh);
    
    for (int sipm = 0; sipm < 4; sipm++) {
      const mu2e::CrvRecoPulse* barPulse = bar->lastPulseInWindow(sipm);
      if (barPulse) { //If we have a valid pulse for the bar
	int color = ((barPulse->GetPulseTime() - ftimeLow) / (ftimeHigh - ftimeLow) * 999);
	bar->SetFillColor(colorPalette[color], sipm);
      }
      else { 
	// Make the SiPM white since no pulses fall within the window
	bar->SetFillColor(kWhite,sipm);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// draw crv
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintCrv(Option_t* Option) {
  int nbars = NBars();
  for (int i = 0; i < nbars; i++) {
    TEvdCrvBar* bar = EvdBar(i);
    bar->Paint(Option);
  }
}

//-----------------------------------------------------------------------------
void TCrvVisNode::PaintXY(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintRZ(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintTZ(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintPhiZ(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintCal(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintVST(Option_t* Option) {
}
//-----------------------------------------------------------------------------
void TCrvVisNode::PaintVRZ(Option_t* Option) {
}

//-----------------------------------------------------------------------------
Int_t TCrvVisNode::DistancetoPrimitive(Int_t px, Int_t py) {
  int              min_dist(9999);
  static TVector3  global;
  
  global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);
  
  //  printf("px,py,X,Y = %5i %5i %10.3f %10.3f\n",px,py,global.X(),global.Y());

  min_dist = DistancetoPrimitiveXY(px, py);
  
  return min_dist;
}

//-----------------------------------------------------------------------------
Int_t TCrvVisNode::DistancetoPrimitiveXY(Int_t px, Int_t py) {

  Int_t min_dist = 9999;
  
  static TVector3 global;
  static TVector3 locrv;
  double          dist;
  
  fClosestObject = NULL;
	
  global.SetXYZ(gPad->AbsPixeltoX(px), gPad->AbsPixeltoY(py), 0);

  double gx = global.X();
  double gy = global.Y();
  
  int          nbar;
  TEvdCrvBar  *evd_bar;
  
  nbar = NBars();
  double gx0, gy0;
  
  printf("gx: %f , gy %f \n", gx, gy);
  for (int ibar = 0; ibar<nbar; ibar++) {
    evd_bar = EvdBar(ibar);
    
    //px0 = gPad->XtoAbsPixel(evd_bar->X0());
    //py0 = gPad->YtoAbsPixel(evd_bar->Y0());
    gx0 = evd_bar->X0();
    gy0 = evd_bar->Y0();
    
    dist = sqrt(pow((gx0 - gx),2) + pow((gy0 - gy),2)); //actual distance
    if (dist < min_dist) {
      min_dist = dist;
      fClosestObject = evd_bar;
      printf("Found bar: %i (%f, %f) to be closest at (%f, %f) dist= %f \n", evd_bar->Bar().asInt(), gx0, gy0, gx, gy, dist);
    }
  }
  return min_dist;
}

//-----------------------------------------------------------------------------
Int_t TCrvVisNode::DistancetoPrimitiveRZ(Int_t px, Int_t py) {
  return 9999;
}

//-----------------------------------------------------------------------------
void TCrvVisNode::Clear(Option_t* Opt) {
  int nbars = NBars();
  for (int ibar = 0; ibar<nbars; ibar++)  {
    EvdBar(ibar)->Clear();
  }
}

//-----------------------------------------------------------------------------
void TCrvVisNode::Print(Option_t* Opt) const {
  TCrvVisNode* node = (TCrvVisNode*) this;

  int nbars = NBars();
  for (int ibar = 0; ibar<nbars; ibar++) {
    TEvdCrvBar* bar = node->EvdBar(ibar);
    bar->Print(Opt);
  }
}

//-----------------------------------------------------------------------------
TEvdCrvBar*  TCrvVisNode::EvdBarWithIndex(int barIndex) {
  TEvdCrvBar* tempBarPtr = 0;

  int nbars = NBars();
  for (int i = 0; i<nbars; i++) {
    tempBarPtr = EvdBar(i);
    if (tempBarPtr->Bar().asInt() == barIndex) {
      return tempBarPtr;
    }
  }
  
  return tempBarPtr; //If it makes it through the loop without returning
}

//-----------------------------------------------------------------------------
void TCrvVisNode::SetTimeWindow(float timeLow, float timeHigh) {
  ftimeLow = timeLow;
  ftimeHigh = timeHigh; 
}

//-----------------------------------------------------------------------------
int TCrvVisNode::getCRVSection(int shieldNumber) {
  int CRVSection = -1;
  switch (shieldNumber) {
  case 0:
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:  CRVSection = 0; break;  //R
  case 6:
  case 7:
  case 8:  CRVSection = 1; break;  //L
  case 9:  CRVSection = 8; break;  //TS T
  case 10:
  case 11:
  case 12: CRVSection = 2; break;  //T
  case 13: CRVSection = 3; break;  //D
  case 14: CRVSection = 4; break;  //U
  case 15: CRVSection = 5; break;  //CU
  case 16: CRVSection = 6; break;  //CD
  case 17: CRVSection = 7; break;  //CT
  case 18:
  case 19:
  case 20: CRVSection = 3; break; //D lower shields - added 05/04/2015
  }
  return CRVSection;
}
