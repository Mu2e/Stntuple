///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCaloRecoDigiBlock__
#define __InitCaloRecoDigiBlock__

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCaloRecoDigiBlock.hh"

namespace stntuple {
class InitCaloRecoDigiBlock : public TStnInitDataBlock {
public:
  art::InputTag   fCaloRecoDigiCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetCaloRecoDigiCollTag (art::InputTag& Tag) { fCaloRecoDigiCollTag = Tag; }
  
  virtual int InitDataBlock    (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks     (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};
}
#endif
