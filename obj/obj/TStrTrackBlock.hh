#ifndef TStrTrackBlock_hh
#define TStrTrackBlock_hh

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TStrTrack.hh"
#include "TBuffer.h"

namespace stntuple {
  class InitStrTrackBlock;
}

class TStrTrackBlock : public TStnDataBlock {
  friend class stntuple::InitStrTrackBlock;
public:
                                        // this is version v1
  int            fNTracks;
  TClonesArray*  fListOfTracks;	//
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TStrTrackBlock();
  virtual ~TStrTrackBlock();
					// ****** accessors

  Int_t       NTracks () { return fNTracks; }

  TStrTrack*  Track(int I) { 
    return (TStrTrack*) fListOfTracks->UncheckedAt(I);
  }

  TClonesArray* GetListOfTracks () { return fListOfTracks ; } 
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TStrTrack*  NewTrack(int ID) {
    return new ((*fListOfTracks)[fNTracks++]) TStrTrack(ID);
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
//  void        ReadV1(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void        Clear(Option_t* opt="") override;
  void        Print(Option_t* opt="") const override;

  ClassDefOverride(TStrTrackBlock,1)
};

#endif
