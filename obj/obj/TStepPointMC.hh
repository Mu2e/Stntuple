///////////////////////////////////////////////////////////////////////////////
// Mu2e Sim Particle - store only primary ones
///////////////////////////////////////////////////////////////////////////////
#ifndef STNTUPLE_TStepPointMC
#define STNTUPLE_TStepPointMC

#include "TMath.h"
#include "TObject.h"
#include "TParticlePDG.h"
#include "TVector3.h"
#include "TBuffer.h"

class TStepPointMC : public TObject {
public:
  int             fVolumeID;
  int             fGenIndex;
  int             fSimID;
  int             fPDGCode;
  int             fCreationCode;
  int             fParentSimID;
  int             fParentPDGCode;
  int             fEndProcessCode;

  float           fEDepTot;
  float           fEDepNio;
  float           fTime;
  float           fStepLength;
  float           fProperTime;          // added in V2;   T(proper)/tau

  TVector3        fPos;			// starting point of the step
  TVector3        fMom;

public:
//------------------------------------------------------------------------------
//  functions
//------------------------------------------------------------------------------
				// ****** constructors and destructor
  TStepPointMC();

  TStepPointMC(int VolumeID    , int GenIndex      ,  
	       int SimID       , int PDGCode       , 
	       int ParentSimID , int ParentPDGCode , 
	       int CreationCode, int EndProcessCode,
	       float EDepTot   , float EDepNio     , 
	       float Time      , float ProperTime  , 
	       float StepLength,
	       float X, float Y, float Z, float Px, float Py, float Pz   );

  virtual ~TStepPointMC();
//-----------------------------------------------------------------------------
// init methods
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int    VolumeID      () const { return fVolumeID;      }
  int    GenIndex      () const { return fGenIndex;      }
  int    SimID         () const { return fSimID;         }
  int    PDGCode       () const { return fPDGCode;       }
  int    ParentSimID   () const { return fParentSimID;   }
  int    ParentPDGCode () const { return fParentPDGCode; }
  int    CreationCode  () const { return fCreationCode;  }
  int    EndProcessCode() const { return fEndProcessCode;}
  float  EDepTot       () const { return fEDepTot;       }
  float  EDepNio       () const { return fEDepNio;       }
  float  Time          () const { return fTime;          }
  float  ProperTime    () const { return fProperTime;    }
  float  SurvivalProb  () const { return exp(-fProperTime); }
  float  StepLength    () const { return fStepLength;       }

  const TVector3*  Pos       () const { return &fPos; }
  const TVector3*  Mom       () const { return &fMom; }
//-----------------------------------------------------------------------------
// overloaded methods of TObject
//-----------------------------------------------------------------------------
  void     Print(Option_t* opt = "") const;
//-----------------------------------------------------------------------------
// schema evolution
//-----------------------------------------------------------------------------
  void     ReadV1(TBuffer &R__b);

  ClassDefOverride(TStepPointMC,2)
};

#endif


