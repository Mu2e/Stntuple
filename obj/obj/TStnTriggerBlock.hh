#ifndef STNTUPLE_TStnTriggerBlock
#define STNTUPLE_TStnTriggerBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE trigger block for Mu2e
//  Author:    P.Murat 
//-----------------------------------------------------------------------------
#include "TClonesArray.h"
#include "Stntuple/obj/TStnDataBlock.hh"
#include "Stntuple/base/TBitset.hh"

class TStnTriggerTable;

class TStnTriggerBlock: public TStnDataBlock {
  friend class StntupleInitTriggerBlock;

public:
//------------------------------------------------------------------------------
//  data members: TSID - fL3Info saved by one copy 
//  I suspect that in versions <=4 last several words of fL3Info are not 
//  restored properly... but the rest - should...
//------------------------------------------------------------------------------
  Short_t   fNPaths;		// # of paths
  Short_t   fNModSeq;		// # of module sequences
  Short_t   fNModules;		// ***** not implemented yet
  Short_t   fExecVersion;
  Short_t   fPhysTable;
  Short_t   fCalibJobSet;
  TBitset   fPaths;		// each bit is 1 is path passed
//------------------------------------------------------------------------------
//  functions: constructors and destructor
//------------------------------------------------------------------------------
public:

  TStnTriggerBlock();
  virtual ~TStnTriggerBlock();
//-----------------------------------------------------------------------------
//  accessors
//-----------------------------------------------------------------------------
					// ****** init methods
  Int_t    Init (Int_t NPaths);
					// ****** accessors

  Int_t    NPaths   () { return fNPaths;   }
  Int_t    NModSeq  () { return fNModSeq;  }
  Int_t    NModules () { return fNModules; }

  Int_t    ExecVersion () { return fExecVersion; }
  Int_t    PhysTable   () { return fPhysTable;   }
  Int_t    CalibJobSet () { return fCalibJobSet; }

  TBitset* Paths       () { return &fPaths; }

  Int_t    PathPassed(Int_t I) { return (fPaths.GetNBits() > I) ? fPaths.GetBit(I) : 0 ; }

//-----------------------------------------------------------------------------
// this could be precalculated
//-----------------------------------------------------------------------------
  Int_t    NPassedPaths();
					// ****** modifiers

  void     SetExecVersion(Int_t Version) { fExecVersion = Version; }
  void     SetNModSeq    (Int_t NModSeq) { fNModSeq     = NModSeq; }
  void     SetPhysTable  (Int_t Table  ) { fPhysTable   = Table  ; }
  void     SetCalibJobSet(Int_t JobSet ) { fCalibJobSet = JobSet ; }

					// for a given trigger 'Table', fills 
					// 'List' with the triggers of level
					// 'Level' (1,2,3) which are fired 
					// for given event and contain pattern
					// 'Name' in their names 
					// function value on return: number of 
					// found triggers

  int      GetListOfPassedTriggers(const TStnTriggerTable* Table, 
				   const char*             Name, 
				   Int_t                   Level, 
				   TObjArray*              List);
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt="") override;
  void Print(Option_t* opt="") const override;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  //  void  ReadV1(TBuffer &R__b);

  ClassDefOverride(TStnTriggerBlock,1)

};

#endif
