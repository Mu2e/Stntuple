#ifndef STNTUPLE_TStrawHitBlock
#define STNTUPLE_TStrawHitBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"

#include "Stntuple/obj/TStrawHit.hh"
#include "Stntuple/obj/TStrWaveform.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

class TStrawHitBlock: public TStnDataBlock {
  friend class stntuple::InitStrawHitBlock;

public:
  int            fNHits;                // number of hits in the straw tracker
  int            fNWaveforms;           // added in V2
  TClonesArray*  fListOfHits;           // list of hits
  TClonesArray*  fListOfWaveforms;      // added in V2, list of waveforms,
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
                                        // constructors and destructor
  TStrawHitBlock();
  virtual ~TStrawHitBlock();
                                        // accessors

  int           NHits     () { return fNHits     ; }
  int           NWaveforms() { return fNWaveforms; }

  TStrawHit*    Hit      (int I) { return (TStrawHit*   ) fListOfHits->UncheckedAt     (I); }
  TStrWaveform* Waveform (int I) { return (TStrWaveform*) fListOfWaveforms->UncheckedAt(I); }
  
  TClonesArray* GetListOfHits     () { return fListOfHits     ; }
  TClonesArray* GetListOfWaveforms() { return fListOfWaveforms; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
                                        //Create hit, increse number of hits

  TStrawHit*    NewHit     (int I) { return new ((*fListOfHits)[fNHits++])           TStrawHit   (I); } 
  TStrWaveform* NewWaveform(int I) { return new ((*fListOfWaveforms)[fNWaveforms++]) TStrWaveform(I); } 
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="") override;
  virtual void Print(Option_t* opt="") const override;

  ClassDefOverride(TStrawHitBlock,2)    // straw hit data block
};

#endif
