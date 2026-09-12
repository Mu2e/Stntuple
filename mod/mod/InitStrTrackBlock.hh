///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitStrTrackBlock__
#define __InitStrTrackBlock__

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TStrTrackBlock.hh"

namespace stntuple {
class InitStrTrackBlock : public TStnInitDataBlock {
public:
  art::InputTag   fTrackCollTag;
  art::InputTag   fTcCollTag;           // time cluster collection
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetTrackCollTag (art::InputTag& Tag) { fTrackCollTag = Tag; }
  void   SetTcCollTag    (art::InputTag& Tag) { fTcCollTag    = Tag; }
  
  virtual int InitDataBlock    (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks     (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};
}
#endif
