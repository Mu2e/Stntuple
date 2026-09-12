#ifndef STNTUPLE_TStepPointMCBlock
#define STNTUPLE_TStepPointMCBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE GENP (generator-level) particle block
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 25 2000
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale 
//  TStepPointMCBlock will be written in split mode, however the TClonesArray 
//  itself may not be split
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "TArrayI.h"
#include "TClonesArray.h"

#include "Stntuple/base/TStnArrayI.hh"
#include "Stntuple/obj/TStepPointMC.hh"
#include "Stntuple/obj/TStnDataBlock.hh"

#include "TSimParticle.hh"
#include "TBuffer.h"

class TStnEvent;

class TStepPointMCBlock: public TStnDataBlock {
  friend class StntupleInitStepPointMCBlock;
protected:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNStepPoints;		// total # of StepPointMC's
  int            fG4Status;		// added in V2 (from MCDataProducts/inc/G4Status.hh)
  int            fNG4Tracks;		// added in V2
  int            fNOverflowSimP;	// added in V2
  int            fNKilledStepLim;	// added in V2
  int            fNKilledFieldProp;	// added in V2
  float          fG4CpuTime;		// added in V2
  float          fG4RealTime;		// added in V2
  TClonesArray*  fListOfStepPoints;
//-----------------------------------------------------------------------------
// transients - parameters (temp solution)
//-----------------------------------------------------------------------------
  int            fGenProcessID;         //! don't save, generated process ID
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TStepPointMCBlock();
  virtual ~TStepPointMCBlock();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int             GenProcessID    () { return fGenProcessID    ; }
  int             NStepPoints     () { return fNStepPoints     ; }
  int             G4Status        () { return fG4Status        ; }
  int             NG4Tracks       () { return fNG4Tracks       ; }
  int             NOverflowSimP   () { return fNOverflowSimP   ; }
  int             NKilledStepLim  () { return fNKilledStepLim  ; }
  int             NKilledFieldProp() { return fNKilledFieldProp; }
  float           G4CpuTime       () { return fG4CpuTime       ; }
  float           G4RealTime      () { return fG4RealTime      ; }

					// `i'-th StepPoint in the list
  TStepPointMC*   StepPointMC(int i) { 
    return (TStepPointMC*) fListOfStepPoints->UncheckedAt(i); 
  }
//-----------------------------------------------------------------------------
//  modifiers
//-----------------------------------------------------------------------------
  TStepPointMC*  NewStepPointMC(int VolumeID, int GenIndex, 
				int SimID   , int PDGCode, 
				int ParentSimID, int ParentPDGCode, 
				int CreationCode, int EndProcessCode,
				float EDepTot, float EDepNio, 
				float Time, float ProperTime, float StepLength,
				float    X, float      Y, float          Z, 
				float   Px, float     Py, float         Pz);

  void           SetGenProcessID(int ID) { fGenProcessID = ID; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "") override;
  virtual void Print(Option_t* opt = "") const override;
//-----------------------------------------------------------------------------
//  I/O and schema evolution
//-----------------------------------------------------------------------------
  void   ReadV1(TBuffer& R__b);

  ClassDefOverride(TStepPointMCBlock,2)	 // TStepPointMC block: Mu2e G4 output
};

#endif


