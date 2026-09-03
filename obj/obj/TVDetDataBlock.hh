#ifndef STNTUPLE_TVDetDataBlock
#define STNTUPLE_TVDetDataBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TVDetHitData.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "TBuffer.h"

class TVDetDataBlock: public TStnDataBlock {
  friend Int_t StntupleInitMu2eVDetDataBlock(TStnDataBlock* block, 
					     AbsEvent*      event, 
					     int            mode);
public:
  Int_t          fNHits;	        // number of hits in the virtual detectors
  TClonesArray*  fListOfHits;		// list of hits
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TVDetDataBlock();
  virtual ~TVDetDataBlock();
					// ****** accessors
  Int_t          NHits    () { return fNHits; }
  TVDetHitData* Hit (int i) { return (TVDetHitData*) fListOfHits->UncheckedAt(i); }
  
  TClonesArray* GetListOfHits () { return fListOfHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        //Create hit, increse number of hits

  TVDetHitData* NewHit() { return new ((*fListOfHits)[fNHits++]) TVDetHitData(); } 
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDefOverride(TVDetDataBlock,1)	// virtual data block
};


#endif


