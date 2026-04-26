#ifndef STNTUPLE_TCrvPulseBlock
#define STNTUPLE_TCrvPulseBlock

#include "TClonesArray.h"

#include "Stntuple/obj/TStnLinkBlock.hh"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/obj/TCrvRecoPulse.hh"
// #include "Stntuple/obj/TCrvCoincidence.hh"

#include "Stntuple/mod/InitStntupleDataBlocks.hh"

class TCrvPulseBlock: public TStnDataBlock {
  friend class StntupleInitCrvPulseBlock;

public:
  int            fNPulses;	             // number of reconstructed pulses
  TClonesArray*  fListOfPulses;	             // list of pulses

  //  int            fNCoincidences;             // number of reconstructed CrvCoincidence's
  //  TClonesArray*  fListOfCoincidences;	     // list of those
  //  TStnLinkBlock* fCoincidencePulseLinks;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TCrvPulseBlock();
  virtual ~TCrvPulseBlock();
					// ****** accessors

  Int_t                   NPulses         () { return fNPulses; }
  TCrvRecoPulse*          Pulse      (int i) { return (TCrvRecoPulse*) fListOfPulses->UncheckedAt(i); }
  TClonesArray*           GetListOfPulses () { return fListOfPulses; }

  //  Int_t                   NCoincidences   () { return fNCoincidences; }
  //  TCrvCoincidence*        Coincidence(int i) { return (TCrvCoincidence*) fListOfCoincidences->UncheckedAt(i); }
  //  TClonesArray*           GetListOfCoincidences() { return fListOfCoincidences; }

  // TStnLinkBlock*          CoincidenceToPulseLinks() { return fCoincidencePulseLinks; }

  // int     NCoincidencePulses   (int Ic        ) const { return fCoincidencePulseLinks->NLinks(Ic); }
  // int     CoincidencePulseIndex(int Ic, int Ip) const { return fCoincidencePulseLinks->Index(Ic,Ip); }
//-----------------------------------------------------------------------------
// modifiers
//-----------------------------------------------------------------------------
  TCrvRecoPulse*          NewPulse() { 
    return new ((*fListOfPulses)[fNPulses++]) TCrvRecoPulse();
  }

  // TCrvCoincidence*        NewCoincidence() { 
  //   return new ((*fListOfCoincidences)[fNCoincidences++]) TCrvCoincidence(); 
  // }

//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="");
  void Print(Option_t* opt="") const;

  ClassDef(TCrvPulseBlock,1)	// CRV reco block
};


#endif
