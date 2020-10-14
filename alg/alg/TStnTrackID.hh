#ifndef murat_TStnTrackID
#define murat_TStnTrackID
///////////////////////////////////////////////////////////////////////////////
// class representing Track ID cuts
// Author:    P. Murat, borrowed from CDF code
// Date:      
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TString.h"

class  TBuffer;

class TStnTrack;
class TH1F;

class TStnTrackID: public TNamed {
public:
//-----------------------------------------------------------------------------
//  try 
//-----------------------------------------------------------------------------
  enum { 
    kNActiveBit        = 0x1 <<  0,     // 0x00000001
    kFitConsBit        = 0x1 <<  1,     // 0x00000002
    kT0Bit             = 0x1 <<  2,     // 0x00000004
    kT0ErrBit          = 0x1 <<  3,     // 0x00000008
    kMomErrBit         = 0x1 <<  4,     // 0x00000010
    kTanDipBit         = 0x1 <<  5,     // 0x00000020
    kD0Bit             = 0x1 <<  6,     // 0x00000040
    kRMaxBit           = 0x1 <<  7,     // 0x00000080
    kTrkQualBit        = 0x1 <<  8,     // 0x00000100
    kChi2DofBit        = 0x1 <<  9,     // 0x00000200
    kDNaBit            = 0x1 << 10,     // 0x00000200
    kFNaBit            = 0x1 << 11      // 0x00000400
  };

  enum { 
    kNFreeInts         =  4,
    kNFreeFloats       =  9
  };

  struct Hist_t {
					// tight object cuts
    TH1F*    fNActive     [5];
    TH1F*    fFitCons     [5];
    TH1F*    fChi2Dof     [5];
    TH1F*    fMomErr      [5];
    TH1F*    fT0Err       [5];
    TH1F*    fT0          [5];
    TH1F*    fDNa         [5];
    TH1F*    fFNa         [5];
    TH1F*    fTanDip      [5];
    TH1F*    fD0          [5];
    TH1F*    fRMax        [5];
    TH1F*    fTrkQual     [5];
					// summary histogram
    TH1F*    fFailedBits;
    TH1F*    fPassed;
  };

protected:
  Int_t      fUseMask;
  Int_t      fMinNActive;
  Int_t      fMaxNActive;
  int        fMaxDNa;
  Int_t      fInteger[kNFreeInts];		// for future

  Float_t    fMinFitCons;		// 
  Float_t    fMinT0;
  Float_t    fMaxT0Err;

  Float_t    fMaxMomErr;
  Float_t    fMinTanDip;
  Float_t    fMaxTanDip;

  Float_t    fMinD0;			// D0 - minimal radius, asymmetric cut
  Float_t    fMaxD0;

  Float_t    fMinRMax;	                // RMax - maximal track radius, asymmetric cut
  Float_t    fMaxRMax;

  Float_t    fMinTrkQual;               // min track quality
  float      fMaxChi2Dof;		// max chi2/DOF (as alternative)
  float      fMinFNa;                   // fraction of active hits

  Float_t    fFloat[kNFreeFloats];	// spare words, added in V5

  void*      fEOR;		// ! end of record
//-----------------------------------------------------------------------------
//  methods
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor

  TStnTrackID(const char* Name = "Anonymous");
  virtual ~TStnTrackID();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Float_t MinFitCons   () const { return fMinFitCons;   }
  Float_t MinT0        () const { return fMinT0;        }
  Int_t   MinNActive   () const { return fMinNActive;   }
  Int_t   MaxNActive   () const { return fMaxNActive;   }
  Float_t MaxT0Err     () const { return fMaxT0Err;     }
  Float_t MaxMomErr    () const { return fMaxMomErr;    }
  Float_t MinTanDip    () const { return fMinTanDip;    }
  Float_t MaxTanDip    () const { return fMaxTanDip;    }
  Float_t MinTrkQual   () const { return fMinTrkQual;   }
  int     MaxDNa       () const { return fMaxDNa    ;   }
  float   MinFNa       () const { return fMinFNa    ;   }
  float   MaxChi2Dof   () const { return fMaxChi2Dof;   }
  Int_t   MinD0        () const { return fMinD0;   }
  Int_t   MaxD0        () const { return fMaxD0;   }
  Int_t   MinRMax      () const { return fMinRMax;   }
  Int_t   MaxRMax      () const { return fMaxRMax;   }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  void    SetMinFitCons  (Float_t FitCons) { fMinFitCons    = FitCons; }
  void    SetMaxChi2Dof  (Float_t Chi2Dof) { fMaxChi2Dof    = Chi2Dof; }
  void    SetMaxDNa      (int     DNa    ) { fMaxDNa        = DNa    ; }
  void    SetMinFNa      (int     FNa    ) { fMinFNa        = FNa    ; }
  void    SetMinT0       (Float_t T0     ) { fMinT0         = T0;      }
  void    SetMinNActive  (Int_t   N      ) { fMinNActive    = N;       }
  void    SetMaxNActive  (Int_t   N      ) { fMaxNActive    = N;       }
  void    SetMaxT0Err    (Float_t T0Err  ) { fMaxT0Err      = T0Err;   }
  void    SetMaxMomErr   (Float_t MomErr ) { fMaxMomErr     = MomErr;  }
  void    SetMinD0       (float   D0     ) { fMinD0         = D0    ;  }
  void    SetMaxD0       (float   D0     ) { fMaxD0         = D0    ;  }
  void    SetMinRMax     (float   R      ) { fMinRMax       = R     ;  }
  void    SetMaxRMax     (float   R      ) { fMaxRMax       = R     ;  }

  void    SetMinTanDip   (Float_t TanDip ) { fMinTanDip     = TanDip;  }
  void    SetMaxTanDip   (Float_t TanDip ) { fMaxTanDip     = TanDip;  }

  void    SetMinTrkQual  (Float_t TrkQual) { fMinTrkQual    = TrkQual; }
  void    SetUseMask     (int     Mask)    { fUseMask       = Mask   ; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  virtual Int_t  IDWord     (TStnTrack* Track);
  virtual Int_t  TightIDWord(TStnTrack* Track);
  virtual Int_t  LooseIDWord(TStnTrack* Track);

                                        // Mode=1: tight ID, =2: loose ID

  void FillHistograms(Hist_t* Hist, TStnTrack* Track, Int_t Mode=1);
//-----------------------------------------------------------------------------
//  overloaded methods of TObject
//-----------------------------------------------------------------------------
  void  Print  (Option_t*  Option = "") const ;
//-----------------------------------------------------------------------------
//  schema evolution
//  2009-12-24: current version - V5
//-----------------------------------------------------------------------------
//   void  ReadV1(TBuffer& R__b);
//   void  ReadV2(TBuffer& R__b);
//   void  ReadV3(TBuffer& R__b);
//   void  ReadV4(TBuffer& R__b);

  ClassDef(TStnTrackID,0)

};

#endif
