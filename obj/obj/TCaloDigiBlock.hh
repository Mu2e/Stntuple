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
  TClonesArray*  fListOfDigis;          //
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:

  TCaloDigiBlock();
  virtual ~TCaloDigiBlock();
					// ****** accessors

  int         NDigis  () { return fNDigis; }

  TCaloDigi*  Digi(int I) { 
    return (TCaloDigi*) fListOfDigis->UncheckedAt(I);
  }

  TClonesArray* GetListOfDigis () { return fListOfDigis ; } 
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCaloDigi*  NewDigi(int I) { 
    return new ((*fListOfDigis)[fNDigis++]) TCaloDigi(I);
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
