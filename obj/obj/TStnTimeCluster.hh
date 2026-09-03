//
// Read the track 
//
// $Id:   $
// $Author:$
// $Date: $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_obj_TStnTimeCluster_hh
#define Stntuple_obj_TStnTimeCluster_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TBuffer.h"

//namespace murat {

namespace mu2e {
  class TimeCluster;
}

class TStnTimeCluster : public TObject {

  enum {
    kNFreeIntsV2   = 10,		// V2
    kNFreeFloatsV2 = 10, 		// V2

    kNFreeIntsV3   =  7,		// V3
    kNFreeFloatsV3 =  9 		// V3
  };

public:
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
  int                       fNHits;
  int                       fHelixSeedIndex;
  int                       fNComboHits;
  int                       fSimID;             // added in V3: Sim ID of the matching MC part 
  int                       fPdgID;             // added in V3: PDG ID of the matching MC part 
  int                       fNHitsSimID;        // added in V3: N(hits) produced by fSimID
  int                       fInt[kNFreeIntsV3]; // provision for future I/O expansion
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
  float                     fT0;    
  float                     fT0Err; 

  float                     fPosX;
  float                     fPosY;
  float                     fPosZ;
 
  float                     fClusterTime;   
  float			    fClusterEnergy; 
  float			    fClusterX;
  float			    fClusterY;
  float			    fClusterZ;

  float                     fMcMom;                 // added in V3: MC truth
  float                     fFloat[kNFreeFloatsV3]; // added in V2: provision for future I/O expansion
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::TimeCluster*    fTimeCluster;          //!
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnTimeCluster(int i = -1);
  ~TStnTimeCluster();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     NHits         () { return  fNHits;  }
  int     NComboHits    () { return  fNComboHits;  }
  int     HelixSeedIndex() { return  fHelixSeedIndex; }

  float   T0            () { return  fT0;     }
  float   T0Err         () { return  fT0Err;  }

  float   posX          () { return  fPosX;}
  float   posY          () { return  fPosY;}
  float   posZ          () { return  fPosZ;}

  float   ClusterTime   () { return fClusterTime;  }
  float   ClusterEnergy () { return fClusterEnergy;}
  float   ClusterX      () { return fClusterX;     }
  float   ClusterY      () { return fClusterY;     }
  float   ClusterZ      () { return fClusterZ;     }

//----------------------------------------------------------------------------
// setters
//----------------------------------------------------------------------------
  void    SetHelixSeedIndex(int I) { fHelixSeedIndex = I; }

//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* Opt = "") override ;
  virtual void Print(Option_t* Opt = "") const override ;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);
  void ReadV2(TBuffer& R__b);

  ClassDefOverride(TStnTimeCluster,3)
};

#endif
