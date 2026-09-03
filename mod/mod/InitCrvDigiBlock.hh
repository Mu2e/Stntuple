///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCrvDigiBlock__
#define __InitCrvDigiBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCrvDigiBlock.hh"

class StntupleInitCrvDigiBlock : public TStnInitDataBlock {
public:
  art::InputTag   fCrvDigiCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetCrvDigiCollTag (art::InputTag& Tag) { fCrvDigiCollTag = Tag; }
  //   void   SetStrawDigiMCCollTag (art::InputTag& Tag) { fStrawDigiMCCollTag = Tag; }
  
  virtual int InitDataBlock    (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks     (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
