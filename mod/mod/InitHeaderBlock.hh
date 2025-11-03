///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __Stntuple_InitHeaderBlock_hh__
#define __Stntuple_InitHeaderBlock_hh__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

namespace stntuple {
  class InitHeaderBlock : public TStnInitDataBlock {
  public:
    art::InputTag   fPbiTag;
    art::InputTag   fShCollTag;
    art::InputTag   fChCollTag;
    art::InputTag   fCalHitCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:

    void   SetChCollTag    (art::InputTag& Tag ) { fChCollTag      = Tag; }
    void   SetShCollTag    (art::InputTag& Tag ) { fShCollTag      = Tag; }
    void   SetCalHitCollTag(art::InputTag& Tag ) { fCalHitCollTag  = Tag; }
    void   SetPbiTag       (art::InputTag& Tag ) { fPbiTag         = Tag; }

    virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode);
    // virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  };
}

#endif
