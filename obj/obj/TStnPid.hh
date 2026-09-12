//
// Read the tracks added to the event by the track finding and fitting code.
//
// $Id: TStnPid.hh,v 1.1 2014/06/13 06:14:48 murat Exp $
// $Author: murat $
// $Date: 2014/06/13 06:14:48 $
//
// Contact person Pavel Murat
//
#ifndef Stntuple_obj_TStnPid_hh
#define Stntuple_obj_TStnPid_hh

// storable objects (data products)
// 

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"

class TStnTrack;

namespace mu2e {
  class AvikPIDNewProduct;
}

class TStnPid : public TObject {

  enum {
    kNFreeInts     = 10,		// V1
    kNFreeFloats   = 10			// V1
  };

public:
//-----------------------------------------------------------------------------
// integers
//-----------------------------------------------------------------------------
  int                       fTrkNumber;     // track number
  int                       fNMatched;
  int                       fNMatchedAll;
  int                       fNUsedOsH;      // for fDrdsSsEle
  int                       fNUsedSsH;      // for fDrdsSsEle
  int                       fNUsedOsD;      // for fSumAvikOsEle
  int                       fInt[kNFreeInts];
//-----------------------------------------------------------------------------
// floats
//-----------------------------------------------------------------------------
  float                     fLogDedxProbEle;
  float                     fLogDedxProbMuo;

  float                     fDrdsVadim;
  float                     fDrdsVadimErr;
  float                     fDrdsOs;
  float                     fDrdsOsErr;
  float                     fDrdsSs;
  float                     fDrdsSsErr;
  
  float                     fSumAvik;
  float                     fSq2Avik;
  float                     fSumAvikOs;

  float                     fFloat[kNFreeFloats];
//-----------------------------------------------------------------------------
// transients
//-----------------------------------------------------------------------------
  const mu2e::AvikPIDNewProduct*  fAvikPid;  //!
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  TStnPid(int i = -1);
  ~TStnPid();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int     TrkNumber()    { return fTrkNumber; }

  float   DrdsVadim   () { return fDrdsVadim;    }
  float   DrdsVadimErr() { return fDrdsVadimErr; }

					// the rest accessors - to come
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* Opt = "") override;
  virtual void Print(Option_t* Opt = "") const override;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
//  void ReadV1(TBuffer& R__b);

  ClassDefOverride(TStnPid,1)
};

#endif
