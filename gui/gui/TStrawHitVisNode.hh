///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TStrawHitVisNode_hh
#define TStrawHitVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"

#else
namespace mu2e {
  class StrawHitCollection;
  class StrawHitPositionCollection;
  class StrawHitFlagCollection;
  //  class PtrStepPointMCVectorCollection;
};
#endif

class TStrawHitVisNode: public TVisNode {
public:
  enum {
    kPickHits     = 0,
    kPickTracks   = 1,
    kPickClusters = 2
  };
  
protected:
  //  TObjArray**    fListOfClusters;

  const mu2e::StrawHitCollection**             fStrawHitColl;
  const mu2e::StrawHitPositionCollection**     fStrawHitPosColl;  //
  const mu2e::StrawHitFlagCollection**         fStrawHitFlagColl; //
  const mu2e::TimeClusterCollection**          fTimeClusterColl;  //
  //  const mu2e::PtrStepPointMCVectorCollection** fMcPtrColl; 
 
  TArc*         fArc;

  TObjArray*    fListOfStrawHits;

  const mu2e::TimeCluster*  fTimePeak;

  Int_t         fDisplayBackgroundHits;
  Int_t         fTimeWindow;
  Int_t         fPickMode;
  Int_t         fUseStereoHits;
  double        fEventTime;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TStrawHitVisNode() {}
  TStrawHitVisNode(const char* Name); 

  virtual ~TStrawHitVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::StrawHitCollection* GetStrawHitColl() { 
    return *fStrawHitColl; 
  }

  const mu2e::StrawHitPositionCollection* GetStrawHitPosColl() { 
    return *fStrawHitPosColl;
  }

  const mu2e::StrawHitFlagCollection* GetStrawHitFlagColl() { 
    return *fStrawHitFlagColl;
  }

  // const mu2e::PtrStepPointMCVectorCollection* GetMcPtrColl() { 
  //   return *fMcPtrColl;
  // }

  int DisplayBackgroundHits() { return fDisplayBackgroundHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetStrawHitColl(const mu2e::StrawHitCollection** Coll) { 
    fStrawHitColl = Coll;
  }

  void SetStrawHitPosColl (const mu2e::StrawHitPositionCollection** Coll) { 
    fStrawHitPosColl = Coll;
  }

  void SetStrawHitFlagColl(const mu2e::StrawHitFlagCollection** Coll) { 
    fStrawHitFlagColl = Coll;
  }

  void SetMcPtrColl(const mu2e::PtrStepPointMCVectorCollection** Coll) { 
    fMcPtrColl = Coll;
  }

  void SetTimeClusterColl(const mu2e::TimeClusterCollection** Coll) { 
    fTimeClusterColl = Coll;
  }

  void  SetPickMode   (Int_t Mode) { fPickMode    = Mode; }

  void  SetDisplayBackgroundHits(Int_t Mode) { fDisplayBackgroundHits = Mode; }

  //  virtual void  Draw    (Option_t* option = "");

  int InitEvent();

  virtual void  Paint   (Option_t* option = "");

  virtual void  PaintXY  (Option_t* option = "") override;
  virtual void  PaintRZ  (Option_t* option = "") override;
  virtual void  PaintTZ  (Option_t* option = "") override;
  virtual void  PaintPhiZ(Option_t* option = "") override;
  virtual void  PaintCal (Option_t* option = "") override;
  virtual void  PaintVST (Option_t* option = "") override;
  virtual void  PaintVRZ (Option_t* option = "") override;

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitive  (Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TStrawHitVisNode,0)
};


#endif
