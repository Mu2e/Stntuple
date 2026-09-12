///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCalDataBlock__
#define __InitCalDataBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCalDataBlock.hh"

namespace stntuple {
class InitCalDataBlock : public TStnInitDataBlock {
public:
  art::InputTag                    fCaloHitCollTag;
  const mu2e::CaloHitCollection*   fCaloHitColl;

  int                              fLastRun;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  InitCalDataBlock();
  
  void        SetCaloHitCollTag(art::InputTag& Tag) { fCaloHitCollTag = Tag; }

  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode) override;
  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode) override;
};
}
#endif
