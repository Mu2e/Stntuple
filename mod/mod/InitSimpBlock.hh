///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef __STNTUPLE_InitSimpBlock__
#define __STNTUPLE_InitSimpBlock__

#include <string.h>

#include "canvas/Utilities/InputTag.h"

#include "Stntuple/obj/TStnInitDataBlock.hh"
#include "Stntuple/obj/TSimpBlock.hh"

class StntupleInitSimpBlock : public TStnInitDataBlock {
public:
  art::InputTag   fSimpCollTag;
  art::InputTag   fStrawHitCollTag;
  art::InputTag   fStrawDigiMCCollTag;
  art::InputTag   fVDHitsCollTag;
  float           fMinSimpEnergy;
  float           fMaxZ;
  int             fGenProcessID;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
public:

  void   SetSimpCollTag       (std::string& Tag ) { fSimpCollTag        = art::InputTag(Tag); }
  void   SetStrawHitCollTag   (std::string& Tag ) { fStrawHitCollTag    = art::InputTag(Tag); }
  void   SetStrawDigiMCCollTag(std::string& Tag ) { fStrawDigiMCCollTag = art::InputTag(Tag); }
  void   SetVDHitsCollTag     (std::string& Tag ) { fVDHitsCollTag      = art::InputTag(Tag); }
  void   SetMinSimpEnergy     (double       MinE) { fMinSimpEnergy      = MinE              ; }
  void   SetMaxZ              (double       MaxZ) { fMaxZ               = MaxZ              ; }
  void   SetGenProcessID      (int          ID  ) { fGenProcessID       = ID                ; }

  virtual int InitDataBlock(TStnDataBlock* Block, AbsEvent* Evt, int Mode);
  //  virtual int ResolveLinks (TStnDataBlock* Block, AbsEvent* Evt, int Mode);

};

#endif
