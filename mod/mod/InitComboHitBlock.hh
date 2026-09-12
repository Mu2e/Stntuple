///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitComboHitBlock__
#define __InitComboHitBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TComboHitBlock.hh"

namespace stntuple {
class InitComboHitBlock : public TStnInitDataBlock {
public:
  art::InputTag   fShCollTag;
  art::InputTag   fChCollTag;
  art::InputTag   fTcCollTag;

  const mu2e::StrawHitCollection*    fShc;
  const mu2e::ComboHitCollection*    fChc;
  const mu2e::TimeClusterCollection* fTcc;

  int             fSaveShLinks;
  
  int             fNTc;
  int             fNSh;
  int             fNCh;
  
  float           fMinDt;
  float           fMinEDep;

  int             fLastRun;

  mu2e::ProditionsHandle<mu2e::TrackerPanelMap>  fTpmh;
  const mu2e::TrackerPanelMap*                   fTpm;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  InitComboHitBlock();
  
  void   SetChCollTag          (art::InputTag& Tag) { fChCollTag          = Tag; }
  void   SetShCollTag          (art::InputTag& Tag) { fShCollTag          = Tag; }
  void   SetTcCollTag          (art::InputTag& Tag) { fTcCollTag          = Tag; }
  
  void   SetMinDt   (float Dt  ) { fMinDt   = Dt  ; }
  void   SetMinEDep (float EDep) { fMinEDep = EDep; }

  bool        CloseEnough  (float Time);
  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) override;
  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode) override;
};
}
#endif
