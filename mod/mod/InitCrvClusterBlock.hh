///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __InitCrvClusterBlock__
#define __InitCrvClusterBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TCrvClusterBlock.hh"

class StntupleInitCrvClusterBlock : public TStnInitDataBlock {
public:
  art::InputTag   fCrvRecoPulseCollTag;
  art::InputTag   fCrvCoincidenceClusterCollTag;
  art::InputTag   fCrvCoincidenceClusterMCCollTag;
  int             fStorePulses;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetCrvRecoPulseCollTag           (std::string& Tag) { fCrvRecoPulseCollTag            = art::InputTag(Tag); }
  void   SetCrvCoincidenceClusterCollTag  (std::string& Tag) { fCrvCoincidenceClusterCollTag   = art::InputTag(Tag); }
  void   SetCrvCoincidenceClusterMCCollTag(std::string& Tag) { fCrvCoincidenceClusterMCCollTag = art::InputTag(Tag); }
  void   SetStorePulses                   (int flag        ) { fStorePulses                    = flag              ; }


  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
