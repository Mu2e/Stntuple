///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitTimeClusterBlock__
#define __InitTimeClusterBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TStnTimeClusterBlock.hh"

namespace stntuple {
  
class InitTimeClusterBlock : public TStnInitDataBlock {
public:
  art::InputTag   fTcCollTag;
  art::InputTag   fShCollTag;
  art::InputTag   fChCollTag;
  art::InputTag   fSdmcCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetTcCollTag  (art::InputTag& Tag) { fTcCollTag   = Tag; }
  void   SetChCollTag  (art::InputTag& Tag) { fChCollTag   = Tag; }
  void   SetShCollTag  (art::InputTag& Tag) { fShCollTag   = Tag; }
  void   SetSdmcCollTag(art::InputTag& Tag) { fSdmcCollTag = Tag; }
  
  virtual int InitDataBlock    (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks     (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};
  
}
#endif
