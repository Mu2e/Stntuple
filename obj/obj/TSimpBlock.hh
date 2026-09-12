#ifndef STNTUPLE_TSimpBlock
#define STNTUPLE_TSimpBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE GENP (generator-level) particle block
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 25 2000
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale 
//  TSimpBlock will be written in split mode, however the TClonesArray 
//  itself may not be split
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "TArrayI.h"
#include "TClonesArray.h"

#include "Stntuple/base/TStnArrayI.hh"
#include "TStnDataBlock.hh"
#include "TSimParticle.hh"
#include "TBuffer.h"

class TStnEvent;

class TSimpBlock: public TStnDataBlock {
  //  friend Int_t StntupleInitMu2eGenpBlock  (TStnDataBlock*, TStnEvent*, int);
  friend class StntupleInitSimpBlock;
protected:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNParticles;		// total # of particles
  TClonesArray*  fListOfParticles;
//-----------------------------------------------------------------------------
// transients - parameters (temp solution
//-----------------------------------------------------------------------------
  int            fGenProcessID;         //! don't save, generated process ID
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TSimpBlock();
  virtual ~TSimpBlock();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  Int_t           NParticles  () { return fNParticles  ; }

					// `i'-th particle in the global list
  TSimParticle*   Particle(int i) { 
    return (TSimParticle*) fListOfParticles->UncheckedAt(i); 
  }

  TSimParticle*   FindParticle(int ID);
//-----------------------------------------------------------------------------
//  modifiers
//-----------------------------------------------------------------------------
					// currently: add LAST particle to the
					// LAST interaction

  TSimParticle*  NewParticle(Int_t ID, Int_t ParentID, Int_t PdgCode, 
			     int CreationCode, int TerminationCode,
			     int StartVolumeIndex, int EndVolumeIndex,
			     int GenProcessID);
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  virtual void Clear(Option_t* opt = "") override;
  virtual void Print(Option_t* opt = "") const override;
//-----------------------------------------------------------------------------
//  I/O and schema evolution
//-----------------------------------------------------------------------------
  ClassDefOverride(TSimpBlock,1)		// TSimBlock: simulated particles
};

#endif


