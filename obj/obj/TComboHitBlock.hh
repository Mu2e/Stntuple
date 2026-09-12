#ifndef STNTUPLE_TComboHitBlock
#define STNTUPLE_TComboHitBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TComboHit.hh"

namespace stntuple {
  class InitComboHitBlock;
}

class TComboHitBlock: public TStnDataBlock {
  friend class stntuple::InitComboHitBlock;

public:
  int            fNHits;                // number of hits in the straw tracker
  TClonesArray*  fListOfHits;           // list of hits
  TStnLinkBlock* fListOfShIndices;      // list of straw hit indices
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
                                        // constructors and destructor
  TComboHitBlock();
  virtual ~TComboHitBlock();
                                        // accessors

  int            NHits     () { return fNHits     ; }

  TComboHit*     Hit      (int I) { return (TComboHit*   ) fListOfHits->UncheckedAt     (I); }
  
  TClonesArray*  ListOfHits     () { return fListOfHits     ; }
  TStnLinkBlock* ListOfShIndices() { return fListOfShIndices; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        // Create hit, increment total number of hits

  TComboHit*    NewHit     (int I) { return new ((*fListOfHits)[fNHits++])           TComboHit   (I); } 
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
//  void ReadV1(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="") override;
  virtual void Print(Option_t* opt="") const override;

  ClassDefOverride(TComboHitBlock,1)    // combo hit data block
};

#endif
