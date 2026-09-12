#ifndef TCaloRecoDigiBlock_hh
#define TCaloRecoDigiBlock_hh

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TCaloRecoDigi.hh"
#include "TBuffer.h"

namespace stntuple {
  class InitCaloRecoDigiBlock;
}

class TCaloRecoDigiBlock : public TStnDataBlock {
  friend class stntuple::InitCaloRecoDigiBlock;
public:
                                        // this is version v1
  int            fNDigis;		//
  TClonesArray*  fListOfCaloRecoDigis;	//
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TCaloRecoDigiBlock();
  virtual ~TCaloRecoDigiBlock();
					// ****** accessors

  Int_t         NDigis         () { return fNDigis; }

  TCaloRecoDigi*  CaloRecoDigi(int I) { 
    return (TCaloRecoDigi*) fListOfCaloRecoDigis->UncheckedAt(I);
  }

  TClonesArray* GetListOfCaloRecoDigis () { return fListOfCaloRecoDigis ; } 
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCaloRecoDigi*  NewCaloRecoDigi(int ID) { 
    return new ((*fListOfCaloRecoDigis)[fNDigis++]) TCaloRecoDigi(ID);
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  // void        ReadV1(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void        Clear(Option_t* opt="") override;
  void        Print(Option_t* opt="") const override;

  ClassDefOverride(TCaloRecoDigiBlock,1)
};

#endif
