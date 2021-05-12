///////////////////////////////////////////////////////////////////////////////
// vis node displays one wedge
///////////////////////////////////////////////////////////////////////////////
#ifndef TComboHitVisNode_hh
#define TComboHitVisNode_hh

#include "Gtypes.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TPad.h"
#include "TArc.h"

#include "Stntuple/base/TVisNode.hh"

#ifndef __CINT__
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#else
namespace mu2e {
  class ComboHitCollection;
};
#endif

class TComboHitVisNode: public TVisNode {
public:
  enum {
    kPickHits     = 0
  };
  
protected:
  const mu2e::ComboHitCollection** fHitColl;
 
  TObjArray*    fListOfHits;

public:
//-----------------------------------------------------------------------------
// constructors and destructor
//-----------------------------------------------------------------------------
  TComboHitVisNode();
  TComboHitVisNode(const char* Name); 

  virtual ~TComboHitVisNode();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const mu2e::ComboHitCollection* GetHitColl() { return *fHitColl; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void SetHitColl(const mu2e::ComboHitCollection** Coll) { fHitColl = Coll; }

  //  virtual void  Draw    (Option_t* option = "");

  int InitEvent();

  virtual void  PaintXY (Option_t* option = "");
  // virtual void  PaintRZ (Option_t* option = "");
  virtual void  PaintTZ (Option_t* option = "");

  //  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);

  virtual Int_t DistancetoPrimitiveXY(Int_t px, Int_t py);
  // virtual Int_t DistancetoPrimitiveRZ(Int_t px, Int_t py);
  virtual Int_t DistancetoPrimitiveTZ(Int_t px, Int_t py);

  //  virtual void   Print(const char* Opt = "") const ; // **MENU**

  ClassDef(TComboHitVisNode,0)
};


#endif
