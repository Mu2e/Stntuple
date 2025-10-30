#ifndef STNTUPLE_TStnHeaderBlock
#define STNTUPLE_TStnHeaderBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE event header
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Oct 31 2000
// 
//-----------------------------------------------------------------------------
#include "TStnDataBlock.hh"
#include "TBuffer.h"

class TStnEvent;

class TStnHeaderBlock : public TStnDataBlock {

  friend Int_t StntupleInitMu2eHeaderBlock     (TStnDataBlock*, AbsEvent* , int);
  friend Int_t StntupleInitMu2eHeaderBlockLinks(TStnDataBlock*, AbsEvent* , int);

  enum {
    kNFreeInts_v2   =  5,               // starting from version 2
    kNFreeFloats_v2 =  5
  };
  enum {
    kNFreeInts   =  4,                  // starting from version 3
    kNFreeFloats =  4
  };

public:
  Int_t             fVersion;
  Int_t             fEventNumber;
  Int_t             fRunNumber;
  Int_t             fSectionNumber;	// section number within the run
  Int_t             fMcFlag;		// MC flag, 0 for real data
  Int_t             fGoodRun;		// run flag
  Int_t             fBrCode;		// 
  Int_t             fGoodTrig;
  Int_t             fTrigWord;		// 
  Int_t             fNTracks;           //
  int               fNStrawHits;        // *** added in V2 
  int               fNCaloHits;         // *** added in V2 - total N reco pulses
  int               fNCRVHits;          // *** added in V2 - total N CRV hits
  Int_t             fCpu;               // packed word with processing time
  int               fNComboHits;        // *** added in V3
  int               fInt[kNFreeInts];   // provision for future I/O expansion

  float             fInstLum;		 // instantaneous luminosity
  float             fMeanLum;		 // *** added in V2 : mean luminosity (MC)
  float             fEventWeight;        // *** added in V4 : event weight
  float             fFloat[kNFreeFloats];// provision for future I/O expansion
  TString           fStnVersion;         // STNTUPLE version, like "dev_243_16"
//------------------------------------------------------------------------------
//  transient data, all persistent should go above
//------------------------------------------------------------------------------
  Int_t             fLastNumber;	//! event/run numbers for the last printed event
  Int_t             fLastRunNumber;	//! 
//------------------------------------------------------------------------------
//  function members
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStnHeaderBlock();
  virtual ~TStnHeaderBlock();
					// ****** accessors 

  int    EventNumber  () const { return fEventNumber;   }
  int    RunNumber    () const { return fRunNumber;     }
  int    SectionNumber() const { return fSectionNumber; }
  int    SubrunNumber () const { return fSectionNumber; }    // Mu2e synonim of the CDF "section number"
  int    McFlag       () const { return fMcFlag;        }
  int    NTracks      () const { return fNTracks;       }
  int    NStrawHits   () const { return fNStrawHits;    }
  int    NComboHits   () const { return fNComboHits;    }

  float InstLum       () const { return fInstLum;       }
  float MeanLum       () const { return fMeanLum;       }
  float LumWeight     () const { return fInstLum/fMeanLum; }

  float CpuTime       () const { return ((fCpu>>8)/10.0);   } // in s
  float CpuSpeed      () const { return ((fCpu&0xFF)/5.0); } // in GHz

  const TString& StnVersion () const { return fStnVersion;    }

					// ****** setters/modifiers

//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void   Clear(Option_t* opt = "");
  virtual void   Print(Option_t* opt = "") const;

					// ****** schema evolution
  void   ReadV1(TBuffer& R__b);


  ClassDef(TStnHeaderBlock,4)	       // Mu2e STNTUPLE event header
};

#endif
