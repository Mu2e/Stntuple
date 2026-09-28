///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCaloClusterBlock__
#define __InitCaloClusterBlock__

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TStnClusterBlock.hh"

namespace stntuple {
  
class InitCaloClusterBlock : public TStnInitDataBlock {
public:
  art::InputTag   fCaloClusterCollTag;
  art::InputTag   fCaloClusterMCCollTag;
  art::InputTag   fCaloHitCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetCaloClusterCollTag      (art::InputTag& Tag) { fCaloClusterCollTag   = Tag; }
  void   SetCaloHitCollTag          (art::InputTag& Tag) { fCaloHitCollTag       = Tag; }
  void   SetCaloClusterMCCollTag    (art::InputTag& Tag) { fCaloClusterMCCollTag = Tag; }

  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Event, int Mode);
  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Event, int Mode);

};
  
}
#endif
