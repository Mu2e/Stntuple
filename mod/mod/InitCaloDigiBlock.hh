///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCaloDigiBlock__
#define __InitCaloDigiBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCaloDigiBlock.hh"

class StntupleInitCaloDigiBlock : public TStnInitDataBlock {
public:
  art::InputTag   fCaloDigiCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetCaloDigiCollTag (art::InputTag& Tag) { fCaloDigiCollTag = Tag; }
  
  virtual int InitDataBlock    (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks     (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
