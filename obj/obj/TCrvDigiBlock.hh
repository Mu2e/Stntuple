 #ifndef STNTUPLE_TCrvDigiBlock
#define STNTUPLE_TCrvDigiBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnLinkBlock.hh"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TCrvDigi.hh"

#include "Stntuple/mod/InitCrvDigiBlock.hh"

class TCrvDigiBlock: public TStnDataBlock {
  friend class StntupleInitCrvDigiBlock;

public:
  int            fNDigis;            // number of reconstructed pulses
  TClonesArray*  fListOfDigis;       // list of pulses
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
                                        // ****** constructors and destructor
  TCrvDigiBlock();
  virtual ~TCrvDigiBlock();
                                        // ****** accessors

  Int_t                   NDigis          () { return fNDigis; }
  TCrvDigi*               Digi       (int i) { return (TCrvDigi*) fListOfDigis->UncheckedAt(i); }
  TClonesArray*           GetListOfDigis  () { return fListOfDigis; }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCrvDigi*          NewDigi() { 
    return new ((*fListOfDigis)[fNDigis++]) TCrvDigi();
  }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="")       override;
  virtual void Print(Option_t* opt="") const override;

  ClassDefOverride(TCrvDigiBlock,1)     // CRV reco block
};

#endif
