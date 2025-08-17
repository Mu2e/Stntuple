///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TMcTruthVisNode_hh
#define TMcTruthVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"
#include "TGraph.h"

#include "Stntuple/gui/TStnVisNode.hh"

#ifndef __CINT__
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

#else

namespace mu2e {
  class SimParticlesWithHits;
  class StepPointMCCollection;
  class GenParticleCollection;
};

#endif

namespace stntuple {
class TMcTruthVisNode: public TStnVisNode {
public:
  enum {
    kPickMcParticles = 0
  };
  
protected:
  const mu2e::StepPointMCCollection**          fSteps; 
  const mu2e::GenParticleCollection**          fGenpColl; 
 
  TArc*         fArc;
  TGraph*       fGraph;

  double        fEventTime;
  double        fTimeWindow;
  double        fMinEnergyDep;
  Int_t         fPickMode;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TMcTruthVisNode() {}
  TMcTruthVisNode(const char* Name); 

  virtual ~TMcTruthVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::StepPointMCCollection* GetStepPointMCCollection() { 
    return *fSteps;
  }

  const mu2e::GenParticleCollection* GetGenpColl() { return *fGenpColl; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetStepPointMCCollection(const mu2e::StepPointMCCollection** Coll) { 
    fSteps = Coll;
  }
  void SetGenpColl(const mu2e::GenParticleCollection** Coll) { fGenpColl = Coll; }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }

  int InitEvent();

  virtual void  PaintXY  (Option_t* option = "");
  virtual void  PaintRZ  (Option_t* option = "");
  virtual void  PaintTZ  (Option_t* option = "");
  virtual void  PaintPhiZ(Option_t* option = "");
  virtual void  PaintCal (Option_t* option = "");
  virtual void  PaintVST (Option_t* option = "");
  virtual void  PaintVRZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  //  ClassDef(TMcTruthVisNode,0)
};
}
#endif
