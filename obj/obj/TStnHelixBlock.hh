#ifndef STNTUPLE_TStnHelixBlock
#define STNTUPLE_TStnHelixBlock
//-----------------------------------------------------------------------------
//  definition of the cluster block for MU2E analysis
//  Author:    G. Pezzullo
//  Date:      April 11 2016
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnHelix.hh"
#include "TBuffer.h"

class TStnHelixBlock: public TStnDataBlock {

  friend class StntupleInitHelixBlock; 
  // friend Int_t StntupleInitMu2eHelixBlock     (TStnDataBlock*, AbsEvent* , int);
  // friend Int_t StntupleInitMu2eHelixBlockLinks(TStnDataBlock*, AbsEvent* , int);
public:
//----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNHelices;
  TClonesArray*  fListOfHelices;
//----------------------------------------------------------------------------
//  functions
//----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnHelixBlock();
  virtual ~TStnHelixBlock();

  TStnHelix* NewHelix() {
    TStnHelix* helixSeed = new ((*fListOfHelices)[fNHelices]) TStnHelix(fNHelices);
    fNHelices++;
    return helixSeed;
  }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NHelices     () { return fNHelices;   }
  TClonesArray*   ListOfHelices() { return fListOfHelices; }

  TStnHelix*   Helix(int I) {
    return (TStnHelix*) fListOfHelices->UncheckedAt(I); 
  }

//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDefOverride(TStnHelixBlock,1)
};

#endif
