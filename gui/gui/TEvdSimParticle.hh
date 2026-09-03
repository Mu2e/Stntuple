///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TEvdSimParticle_hh
#define TEvdSimParticle_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TParticlePDG.h"

#include "Stntuple/obj/TSimParticle.hh"

namespace mu2e {
  class StrawGasStep;
}

namespace stntuple {
class TEvdSimParticle: public TObject {
public:
  int                       fNumber;

  TSimParticle*             fSimp;
  const mu2e::StrawGasStep* fS1;	// first step on the trajectory
  const mu2e::StrawGasStep* fS2;	// last  step on the trajectory, 
					// not sure how to deal with turning particles

  TObjArray*                fListOfHits;
  TEllipse*                 fEllipse;
  TLine*                    fLineXY;	// for zero field
  TLine*                    fLineTZ;
  TParticlePDG*             fParticlePDG;
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TEvdSimParticle();
  TEvdSimParticle(int Number, TSimParticle* TSimp, const mu2e::StrawGasStep* S1, const mu2e::StrawGasStep* S2); 

  virtual ~TEvdSimParticle();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TSimParticle*       SimParticle() { return fSimp; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  virtual int DistancetoPrimitive    (int px, int py);
  virtual int DistancetoPrimitiveXY  (int px, int py);
  virtual int DistancetoPrimitiveRZ  (int px, int py);
  virtual int DistancetoPrimitiveTZ  (int px, int py);
  virtual int DistancetoPrimitivePhiZ(int px, int py);

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py)
//-----------------------------------------------------------------------------
// drawing functions
//-----------------------------------------------------------------------------
  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void  Paint(Option_t* Opt = "");
  virtual void  Clear(Option_t* Opt = "");
  virtual void  Print(Option_t* Opt = "") const ; // **MENU**

  ClassDefOverride(stntuple::TEvdSimParticle,0)
};

}

#endif
