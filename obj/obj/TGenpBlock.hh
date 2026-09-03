#ifndef STNTUPLE_TGenpBlock
#define STNTUPLE_TGenpBlock
//-----------------------------------------------------------------------------
//  definition of the STNTUPLE GENP (generator-level) particle block
//  Author:    Pasha Murat (CDF/FNAL)
//  Date:      Nov 25 2000
// 
//  the implementation is a trade-off between the performance in split and 
//  non-split modes. I'm assuming that on a long term time scale 
//  TGenpBlock will be written in split mode, however the TClonesArray 
//  itself may not be split
//  Also one wants to control the streamer of the class, having a wrapper
//  around TClonesArray seems to be a reasonable compromise here
//-----------------------------------------------------------------------------

#include "TArrayI.h"
#include "TClonesArray.h"

#include "Stntuple/base/TStnArrayI.hh"
#include "TStnDataBlock.hh"
#include "TGenParticle.hh"
#include "TBuffer.h"

class TGenpBlock: public TStnDataBlock {
  friend class StntupleInitGenpBlock;
protected:
//-----------------------------------------------------------------------------
//  data members
//-----------------------------------------------------------------------------
  Int_t          fNParticles;		// total # of particles
  int            fGenProcessID;         // ID of the generated signal process
  float          fWeight;               // generation weight (!= 1 for RPC)
  float          fGenEnergy;            // generation energy if not in the event
  TClonesArray*  fListOfParticles;
//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
					// ****** constructors and destructor
  TGenpBlock();
  virtual ~TGenpBlock();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int             NParticles  () { return fNParticles;   }
  int             GenProcessID() { return fGenProcessID; }
  float           Weight      () { return fWeight;       }
  float           GenEnergy   () { return fGenEnergy;    }

					// `i'-th particle in the global list
  TGenParticle*   Particle(int i) { 
    return (TGenParticle*) fListOfParticles->UncheckedAt(i); 
  }

//-----------------------------------------------------------------------------
//  modifiers
//-----------------------------------------------------------------------------
					// currently: add LAST particle to the
					// LAST interaction

  TGenParticle*  NewParticle(Int_t PdgCode, Int_t GenProcessID, 
			     int   m1, Int_t m2, Int_t d1, Int_t d2, 
			     float px, Float_t py, Float_t pz, Float_t e,
			     float vx, Float_t vy, Float_t vz, Float_t t,
			     float ProperTime);

  void           SetGenProcessID(int ID) { fGenProcessID = ID; }
//-----------------------------------------------------------------------------
// overloaded functions of TObject
//-----------------------------------------------------------------------------
  void Clear(Option_t* opt = "");
  void Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
//  I/O and schema evolution
//-----------------------------------------------------------------------------
  void ReadV1(TBuffer& R__b);
  void ReadV2(TBuffer& R__b);

  ClassDefOverride(TGenpBlock,3)		// GENP block: output of MC event generators
};

#endif


