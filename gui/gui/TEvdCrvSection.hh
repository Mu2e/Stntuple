///////////////////////////////////////////////////////////////////////////////
// Vis node - One section of CRV
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdCrvSection_hh
#define TEvdCrvSection_hh

#include "Rtypes.h" // (?)  - Gtypes is obsolete and contains nothing
#include "TNamed.h"
#include "TBox.h"
#include "TGeoBBox.h"

#ifndef __CINT__

#include "Offline/CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"

#else

namespace mu2e {
  class  CRSScintillatorShield;
};

#endif

class TEvdCrvSection : public TNamed {
public:
  //const mu2e::CRSScintillatorShield*  fShield;
  int       fSectionID;
  int       fNModules;
  int       fNLayers;
  int       fNBars;              // per layer
  int       fFirstBar;
  int       fBarIx;	      // orientation of the bar
  int       fBarIy;
  int       fBarIz;
  float     fBarDx;
  float     fBarDy;
  float     fBarDz;
  TBox*	    fBox;
  TGeoBBox* fBarShape;
public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdCrvSection() {}
  TEvdCrvSection(int SectionID, const char* Name = "");

  virtual ~TEvdCrvSection();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
//const mu2e::CRSScintillatorShield*	Shield()	{ return fShield; }

  TGeoBBox* BarShape() { return fBarShape; }
  
  int SectionID()	{ return fSectionID; }

  //int NBars() { return fShield->(); } - No current support for getting the number of bars in a given section easily
  //-----------------------------------------------------------------------------
  // modifiers
  //-----------------------------------------------------------------------------

  //  virtual void  Draw    (Option_t* option = "");

  virtual void  Paint(Option_t* option = "") override;
  //virtual void  PaintXY(Option_t* Option = "");
  //virtual void  PaintRZ(Option_t* Option = "");
  //virtual void  PaintCal(Option_t* Option = "");
  virtual void  PaintCrv(Option_t* Option = "");

  int   InitEvent() { return 0; }

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py) override;
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDefOverride(TEvdCrvSection, 0)
};


#endif
