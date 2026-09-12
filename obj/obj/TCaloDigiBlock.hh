#ifndef TCaloDigiBlock_hh
#define TCaloDigiBlock_hh

#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/mod/InitStntupleDataBlocks.hh"
#include "Stntuple/obj/TCaloDigi.hh"
#include "TBuffer.h"

class TCaloDigiBlock : public TStnDataBlock {
  friend Int_t StntupleInitMu2eCaloDigiBlock(TStnDataBlock*, AbsEvent*, int);
public:
                                        // this is version v1
  int            fNDigis;		//
  TClonesArray*  fListOfCaloDigis;	//
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TCaloDigiBlock();
  virtual ~TCaloDigiBlock();
					// ****** accessors

  Int_t         NDigis         () { return fNDigis; }

  TCaloDigi*  CaloDigi(int I) { 
    return (TCaloDigi*) fListOfCaloDigis->UncheckedAt(I);
  }

  TClonesArray* GetListOfCaloDigis () { return fListOfCaloDigis ; } 
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCaloDigi*  NewCaloDigi() { 
    return new ((*fListOfCaloDigis)[fNDigis++]) TCaloDigi();
  }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void        ReadV1(TBuffer& R__b);
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void        Clear(Option_t* opt="") override;
  void        Print(Option_t* opt="") const override;

  ClassDefOverride(TCaloDigiBlock,1)
};

#endif
