#ifndef STNTUPLE_TStnTimeClusterBlock
#define STNTUPLE_TStnTimeClusterBlock
//-----------------------------------------------------------------------------
//  definition of the cluster block for MU2E analysis
//  Author:    G. Pezzullo
//  Date:      April 11 2016
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnLinkBlock.hh"
#include "Stntuple/obj/TStnTimeCluster.hh"
#include "TBuffer.h"

namespace stntuple {
  class InitTimeClusterBlock;
}

class TStnTimeClusterBlock: public TStnDataBlock {
  friend class stntuple::InitTimeClusterBlock;
public:
//----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNTimeClusters;
  TClonesArray*  fListOfTimeClusters;
  TStnLinkBlock* fListOfChLinks;        // added in V2
//----------------------------------------------------------------------------
//  functions
//----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnTimeClusterBlock();
  virtual ~TStnTimeClusterBlock();

  TStnTimeCluster* NewTimeCluster() {
    TStnTimeCluster* tc = new ((*fListOfTimeClusters)[fNTimeClusters]) TStnTimeCluster();
    fNTimeClusters++;
    return tc;
  }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t              NTimeClusters     () { return fNTimeClusters;   }
  TClonesArray*      ListOfTimeClusters() { return fListOfTimeClusters; }

  TStnTimeCluster*   TimeCluster(int I) {
    return (TStnTimeCluster*) fListOfTimeClusters->UncheckedAt(I); 
  }

  TStnLinkBlock*     ListOfChLinks() { return fListOfChLinks; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="") override;
  virtual void Print(Option_t* opt="") const override;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer &R__b);
  
  ClassDefOverride(TStnTimeClusterBlock,2)
};

#endif
