#ifndef murat_TStnClusterBlock
#define murat_TStnClusterBlock
//-----------------------------------------------------------------------------
//  definition of the cluster block for MU2E analysis
//  Author:    Pavel Murat (Fermilab)
//  Date:      March 07 2013
//-----------------------------------------------------------------------------
#include "TClonesArray.h"

#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TStnCluster.hh"
#include "Stntuple/obj/TStnLinkBlock.hh"
#include "TBuffer.h"

namespace stntuple {
  class InitCaloClusterBlock;
}

class TStnClusterBlock: public TStnDataBlock {
  friend class stntuple::InitCaloClusterBlock;
public:
//----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNClusters;
  TClonesArray*  fListOfClusters;
  TStnLinkBlock* fListOfHitLinks;             // == added in V2 ==
//----------------------------------------------------------------------------
//  functions
//----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnClusterBlock();
  virtual ~TStnClusterBlock();

  TStnCluster* NewCluster() {
    TStnCluster* cl = new ((*fListOfClusters)[fNClusters]) TStnCluster(fNClusters);
    fNClusters++;
    return cl;
  }
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NClusters     () { return fNClusters;   }
  TClonesArray*   ListOfClusters() { return fListOfClusters; }

  TStnCluster*    Cluster(int I) {
    return (TStnCluster*) fListOfClusters->UncheckedAt(I); 
  }

  TStnLinkBlock*  ListOfHitLinks() { return fListOfHitLinks;     }
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void     ReadV1(TBuffer &R__b);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt="")       override;
  virtual void Print(Option_t* opt="") const override;

  ClassDefOverride(TStnClusterBlock,2)
};

#endif
