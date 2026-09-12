#ifndef STNTUPLE_TTrackStrawHitBlock
#define STNTUPLE_TTrackStrawHitBlock

#include "TClonesArray.h"

#include "Stntuple/base/TStnArrayI.hh"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TTrackStrawHit.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "TBuffer.h"

class TTrackStrawHitBlock: public TStnDataBlock {
  friend class stntuple::InitTrackStrawHitBlock;

public:
  Int_t          fNTracks;	        // total number of tracks
  Int_t          fNStrawHits;
  TStnArrayI*    fNTrackHits;           // number of hits per track 
  TStnArrayI*    fFirst;                // index of the first hit of the i-th track in the common list
  TClonesArray*  fListOfHits;		// list of hits
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TTrackStrawHitBlock();
  virtual ~TTrackStrawHitBlock();
					// ****** accessors
  Int_t          NTracks  ()      { return fNTracks;  }
  Int_t          NTrackHits(int i) { return fNTrackHits->At(i); }
  Int_t          First    (int i) { return fFirst->At(i); } 

  TTrackStrawHit* Hit (int i) { return (TTrackStrawHit*) fListOfHits->UncheckedAt(i); }

  TTrackStrawHit* Hit (int ITrack, int I) { 
    return (TTrackStrawHit*) fListOfHits->UncheckedAt(fFirst->At(ITrack)+I); 
  }
  
  TClonesArray* GetListOfHits () { return fListOfHits; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        // create hit, increse number of hits

  TTrackStrawHit* NewHit() { return new ((*fListOfHits)[fNStrawHits++]) TTrackStrawHit(); } 
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="") override;
  void Print(Option_t* opt="") const override;

  ClassDefOverride(TTrackStrawHitBlock,1)	// track straw hit data block
};


#endif


