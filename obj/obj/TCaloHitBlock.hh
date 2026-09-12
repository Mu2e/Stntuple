#ifndef TCaloHitBlock_hh
#define TCaloHitBlock_hh

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "TCaloHit.hh"
#include "TBuffer.h"
namespace stntuple {
  class InitCaloHitBlock;
}

class TCaloHitBlock : public TStnDataBlock {
  friend class stntuple::InitCaloHitBlock;
  
public:
  int            fNHits;		// number of hit crystals
  float          fEDep[2];              // energy dep per disk
  TClonesArray*  fListOfCaloHits;	// list of crystal hit data 
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TCaloHitBlock();
  virtual ~TCaloHitBlock();
					// ****** accessors

  Int_t         NHits           () { return fNHits; }

  TCaloHit*      Hit(int I) { 
    return (TCaloHit*) fListOfCaloHits->UncheckedAt(I);
  }

  TClonesArray* GetListOfCaloHits () { return fListOfCaloHits ; } 
//-----------------------------------------------------------------------------
// these routines do loops - be careful using them if you need performance
//-----------------------------------------------------------------------------
//   TCaloHitData*  GetTowerByKey(Int_t Key) {
//     return Tower(TCaloHitData::IEta(Key), TCaloHitData::IPhi(Key));
//   }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCaloHit*  NewCaloHit(int ID) { 
    return new ((*fListOfCaloHits)[fNHits++]) TCaloHit(ID);
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  // void        ReadV1(TBuffer& R__b);
  // void        ReadV2(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void        Clear(Option_t* opt="") override;
  virtual void        Print(Option_t* opt="") const override;

  ClassDefOverride(TCaloHitBlock,1)
};

#endif
