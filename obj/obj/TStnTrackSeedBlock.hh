#ifndef murat_TStnTrackSeedBlock
#define murat_TStnTrackSeedBlock
//-----------------------------------------------------------------------------
//  definition of the cluster block for MU2E analysis
//  Author:    G. Pezzullo
//  Date:      April 11 2016
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnTrackSeed.hh"

class TStnTrackSeedBlock: public TStnDataBlock {
  friend class StntupleInitTrackSeedBlock;
public:
//----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNTrackSeeds;
  TClonesArray*  fListOfTrackSeeds;
//----------------------------------------------------------------------------
//  functions
//----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnTrackSeedBlock();
  virtual ~TStnTrackSeedBlock();

  TStnTrackSeed* NewTrackSeed() {
    TStnTrackSeed* trackSeed = new ((*fListOfTrackSeeds)[fNTrackSeeds]) TStnTrackSeed();
    fNTrackSeeds++;
    return trackSeed;
  }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NTrackSeeds     () { return fNTrackSeeds;   }
  TClonesArray*   ListOfTrackSeeds() { return fListOfTrackSeeds; }

  TStnTrackSeed*   TrackSeed(int I) {
    return (TStnTrackSeed*) fListOfTrackSeeds->UncheckedAt(I); 
  }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDefOverride(TStnTrackSeedBlock,1)
};

#endif
