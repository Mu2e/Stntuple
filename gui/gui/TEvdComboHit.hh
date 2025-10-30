///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef Stntuple_gui_TEvdComboHit_hh
#define Stntuple_gui_TEvdComboHit_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TMarker.h"

#include "canvas/Persistency/Provenance/ProductID.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"

namespace mu2e {
  class ComboHit;
  class StrawDigiMC;
  class SimParticle;
}

namespace stntuple {

class TEvdComboHit: public TObject {
public:

  enum { 
    kInTimeBit     = 0x1 << 0,
    kConversionBit = 0x1 << 1
  };

  enum { 
    kStereo      = 0x00000001, // <<  0
    kEnergySel   = 0x00000002, // <<  1
    kRadSel      = 0x00000004, // <<  2
    kTimeSel     = 0x00000008, // <<  3
    kBkgClust    = 0x00000020, // <<  5,
    kBkg         = 0x00000040, // <<  6,
    kIsolated    = 0x00000080, // <<  7,
    kOutlier     = 0x00000100, // <<  8,
    kOther       = 0x00000200, // <<  9,
    kTdiv        = 0x00000400, // << 10,
    kClust       = 0x00000800, // << 11,
    kCalosel     = 0x00001000, // << 12,
    kStrawxtalk  = 0x00002000, // << 13
    kElextalk    = 0x00004000, // << 14
    kTrksel      = 0x00008000, // << 15
    kActive      = 0x00010000, // << 16
    kDoca        = 0x00020000, // << 17
    kResolvedphi = 0x00040000, // << 18
    kCalopresel  = 0x00080000, // << 19
    kIntime      = 0x00100000, // << 20
    kPanelcombo  = 0x00200000, // << 21
    kTrack       = 0x00400000, // << 22
    kDead        = 0x00800000, // << 23
    kNoisy       = 0x01000000, // << 24
  };

protected:
					  // backward pointers to the reconstruction objects
  const mu2e::ComboHit*      fHit;
  const mu2e::SimParticle*   fSim;        // 
  const mu2e::StrawGasStep*  fStep;
  TVector3                   fPos;	  // position in 3D, Z=zwire
  TVector2                   fDir;        // direction of the straw
  double                     fSigW;       // error in the wire direction
  double                     fSigR;       // error in radial direction
  TLine                      fLineW;	  // paint on XY view
  TLine                      fLineR;
  TMarker                    fTZMarker;   // paint on TZ view
  TMarker                    fPhiZMarker; // paint on PhiZ view
  int                        fSimID;
  int                        fPdgID;
  int                        fMotherPdgID;
  float                      fP;
  float                      fPz;
  int                        fColor;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdComboHit();

  TEvdComboHit(const mu2e::ComboHit*     Hit,
	       const mu2e::SimParticle*  Sim,
	       const mu2e::StrawGasStep* Step,
	       int                       MotherPdgID,
	       float                     P,
	       float                     Pz);

  virtual ~TEvdComboHit();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TVector3*              Pos()           { return &fPos; }
  TVector2*              Dir()           { return &fDir; }
  const mu2e::ComboHit*  ComboHit()      { return fHit;  }

  float                  Z()             { return fPos.Z(); }
  float                  T()             { return fHit->time(); }
  float                  correctedTime() { return fHit->correctedTime(); }
  TLine*                 LineW()         { return &fLineW ; }
  TLine*                 LineR()         { return &fLineR ; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  // void SetMask (int Mask ) { fMask = Mask ;}

  void SetColor(int Color) { 
    fLineW.SetLineColor(Color); 
    fLineR.SetLineColor(Color); 
  }

  void SetLineWidth(int W) { 
    fLineW.SetLineWidth(W); 
    fLineR.SetLineWidth(W); 
  }

  virtual void  PaintXY    (Option_t* option = "");
  virtual void  PaintTZ    (Option_t* Option = "");
  // virtual void  Paint      (Option_t* option = "");
  virtual void  PaintPhiZ  (Option_t* option = "");
  // virtual void  PaintCal   (Option_t* option = "");

  // virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  // virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);

  // virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  // virtual Int_t DistancetoPrimitiveTZ(Int_t px, Int_t py);

  virtual void   PrintHeader()                 const ;  // *MENU*
  virtual void   PrintMe()                     const ;  // *MENU*
  virtual void   Print  (const char* Opt = "") const ;  // *MENU*

  ClassDef(stntuple::TEvdComboHit,0)
};

}
#endif
