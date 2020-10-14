///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitStrawDataBlock__
#define __InitStrawDataBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TStrawDataBlock.hh"

class StntupleInitStrawDataBlock : public TStnInitDataBlock {
public:
  art::InputTag   fStrawHitCollTag;
  art::InputTag   fStrawDigiMCCollTag;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetStrawHitCollTag    (std::string& Tag) { fStrawHitCollTag    = art::InputTag(Tag); }
  void   SetStrawDigiMCCollTag (std::string& Tag) { fStrawDigiMCCollTag = art::InputTag(Tag); }
  
  virtual int InitDataBlock  (TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  //  virtual int ResolveLinks   (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
